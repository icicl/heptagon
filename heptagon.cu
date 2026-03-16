#include <cmath>
#include <stdint.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>


#include "hep.h"

double time_us() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ((double)ts.tv_sec * 1000000.0) + ((double)ts.tv_nsec / 1000.0);
}

__device__ uint32_t index(uint32_t a, uint32_t b) { // get the index into the flattened 1D matrix of distances
    return a*(a-1)/2 + b;
}

__device__ Node* make_node(uint32_t left, uint32_t right, float dist, Node* node_data, uint32_t* next_node_data_idx) {
    uint32_t id = *next_node_data_idx;
    node_data[*next_node_data_idx] = (Node){id, left, right, UINT32_MAX, dist};
    if (left != UINT32_MAX) {
        node_data[left].parent = id;
    }
    if (right != UINT32_MAX) {
        node_data[right].parent = id;
    }
    (*next_node_data_idx)++;
    return node_data + id;
}


uint32_t next_compaction_n(uint32_t N) {
    return N < 400 ? 0 : N*7/8;
}

#define COMPACT_NT 1024
__global__ void nj_compact(uint32_t N, float* D, float* D_tmp, float* r, Node** nodes) {
    uint32_t ii, i_sparse, i_dense, idx_sparse, idx_dense, j_sparse, j_dense, nxt_i_dense;
    uint32_t tx = threadIdx.x;
    idx_sparse = 0;
    idx_dense = 0;
    i_dense = 0;
    nxt_i_dense = 0;
    __shared__ uint32_t idx_dense_s[COMPACT_NT];

    for (ii=0; ii < N; ii += COMPACT_NT) {
        i_sparse = ii + tx;
        if (tx == 0) {
            i_dense = nxt_i_dense;
            idx_dense = i_dense * (i_dense - 1) / 2;
            for (uint32_t i=i_sparse; i<i_sparse+COMPACT_NT && i<N; i++) {
                idx_dense_s[i - i_sparse] = idx_dense;
                if (r[i] == -INFINITY) continue;
                idx_dense += i_dense;
                i_dense++;
            }
            nxt_i_dense = i_dense;
        }
        __syncthreads();
        if (i_sparse >= N || r[i_sparse] == -INFINITY) continue;
        idx_dense = idx_dense_s[tx];
        idx_sparse = i_sparse * (i_sparse - 1) / 2;
        j_dense = 0;
        for (j_sparse=0; j_sparse<i_sparse; j_sparse++) {
            if (r[j_sparse] == -INFINITY) continue;
            D_tmp[idx_dense + j_dense] = D[idx_sparse + j_sparse];
            j_dense++;
        }
    }
    __syncthreads();
    if (tx == 0) {
        i_dense = 0;
        for (i_sparse=0; i_sparse < N; i_sparse++) {
            if (r[i_sparse] == -INFINITY) continue;
            r[i_dense] = r[i_sparse];
            nodes[i_dense] = nodes[i_sparse];
            i_dense++;
        }
    }
}


__global__ void nj_setup(uint32_t N_full, Node* node_data_out, Node** nodes, uint32_t* next_node_data_idx) {
    *next_node_data_idx = 0;
    for (uint32_t i=0; i<N_full; i++) nodes[i] = make_node(UINT32_MAX, UINT32_MAX, 0, node_data_out, next_node_data_idx);
}

struct min_data {
    double dist;
    uint32_t i;
    uint32_t j;
};
typedef struct min_data min_data;

#define ARGMIN_NT 256
#define ARGMIN_NB 256
__global__ void nj_argmin(uint32_t N, uint32_t n, const float* D, const float* r, min_data* best_ij) {
    uint32_t tx = threadIdx.x, 
             bx = blockIdx.x;

    uint32_t i, j, best_i, best_j, s;
    float best, dist;

    __shared__ min_data mins[ARGMIN_NT/32];

    best = FLT_MAX;
    for (i=bx; i<N; i+=ARGMIN_NB) {
        if (r[i] == -INFINITY) continue;
        for (j=tx; j<i; j+=ARGMIN_NT) {
            if (r[j] == -INFINITY) continue;
            dist = D[index(i,j)] - r[i] - r[j];
            if (dist < best) {
                best = dist;
                best_i = i;
                best_j = j;
            }
        }
    }

    for (s=16; s>0; s>>=1) { // reduce w/i warp
        dist = __shfl_down_sync(0xFFFFFFFF, best, s);
        i = __shfl_down_sync(0xFFFFFFFF, best_i, s);
        j = __shfl_down_sync(0xFFFFFFFF, best_j, s);
        if (dist < best) {
            best = dist;
            best_i = i;
            best_j = j;
        }
    }
    if (tx%32 == 0) {
        mins[tx/32] = (min_data){best, best_i, best_j};
    }
    __syncthreads();
    for (s = (ARGMIN_NT/2)/32; s>0; s>>=1) {
        if (tx < s) {
            if (mins[tx+s].dist < mins[tx].dist) {
                mins[tx] = mins[tx+s];
            }
        }
        __syncthreads();
    }
    best_ij[bx] = mins[0];
}

#define UPDATE_NT 512
__global__ void nj_update(uint32_t N, uint32_t n, float* D, float* r, Node* node_data_out, Node** nodes, uint32_t* next_node_data_idx, min_data* best_ij) {
    uint32_t tx = threadIdx.x, 
             bx = blockIdx.x;
    uint32_t best_i, best_j, m, idx;
    float dist_i, dist_j, dij, dmi, dmj;
    __shared__ float rjs[UPDATE_NT];

    for (m=ARGMIN_NB/2; m>0; m>>=1) {
        if (tx < m) {
            if (best_ij[tx+m].dist < best_ij[tx].dist) {
                best_ij[tx] = best_ij[tx+m];
            }
        }
        __syncthreads();
    }
    best_i = best_ij[0].i;
    best_j = best_ij[0].j;

    dij = D[index(best_i, best_j)];
    dist_i = (dij + r[best_i] - r[best_j]) * 0.5;
    dist_j = dij - dist_i;
    if (dist_i < 0) dist_i = 0; // negative branch lengths not allowed
    if (dist_j < 0) dist_j = 0; // negative branch lengths not allowed

    if (tx == 0) {
        nodes[best_i]->dist = dist_i;
        nodes[best_j]->dist = dist_j;
        nodes[best_i] = make_node(nodes[best_i]->id, nodes[best_j]->id, 0, node_data_out, next_node_data_idx);
        nodes[best_j] = NULL;
        r[best_j] = -INFINITY;
    }
    __syncthreads();

    float rj = 0.0;
    for (m=tx; m<N; m+=UPDATE_NT) {
        if (nodes[m] == 0 || m == best_i) continue;
        dmj = (m > best_j) ? D[index(m, best_j)] : D[index(best_j, m)];
        idx = (m > best_i) ? index(m, best_i) : index(best_i, m);
        dmi = D[idx];
        D[idx] = (dmi + dmj - dij) * 0.5;
        r[m] = ((r[m] * (n - 2.0)) - dmi - dmj + D[idx]) / (n - 3.0); // subtract distances to i and j, add distance to new merged node, and renormalize
        rj += D[idx];
    }
    rjs[tx] = rj;
    __syncthreads();
    for (m=UPDATE_NT/2; m>0; m>>=1) {
        if (tx < m) rjs[tx] += rjs[tx+m];
        __syncthreads();
    }
    if (tx == 0) r[best_i] = rjs[0] / (n - 3.0);
}

#define ADD_NT 256
__global__ void nj_add(uint32_t n, float* D_row, float* r, Node* node_data_out, Node** nodes, uint32_t* next_node_data_idx) {
    uint32_t best_i, o_best_i, i, j;
    const uint32_t tx = threadIdx.x;
    float best, o_best, q, rn;
    __shared__ float s_best[ADD_NT/32];
    __shared__ uint32_t s_best_i[ADD_NT/32];
    best = FLT_MAX;
    for (i=tx; i<n; i+=ADD_NT) {
        q = D_row[i];
        if (q < best) {
            best = q;
            best_i = i;
        }
    }
    for (i=16; i>0; i>>=1) {
        o_best = __shfl_down_sync(0xFFFFFFFF, best, i);
        o_best_i = __shfl_down_sync(0xFFFFFFFF, best_i, i);
        if (o_best < best) {
            best = o_best;
            best_i = o_best_i;
        }
    }
    if (tx % 32 == 0) {
        s_best[tx/32] = best;
        s_best_i[tx/32] = best_i;
    }
    __syncthreads();
    for (i=(ADD_NT/32)/2; i>0; i>>=1) {
        if (tx < i) {
            if (s_best[tx+i] < s_best[tx]) {
                s_best[tx] = s_best[tx+i];
                s_best_i[tx] = s_best_i[tx+i];
            }
        }
        __syncthreads();
    }
    best_i = s_best_i[0];
    best = D_row[best_i];

    i = node_data_out[best_i].parent;
    if (i == UINT32_MAX) { // if a top-level node
        for (j=tx; j<n; j+=ADD_NT) {
            if (nodes[j] != NULL && nodes[j]->id == best_i) {
                nodes[j] = make_node(best_i, n, best*0.5, node_data_out, next_node_data_idx);
            }
        }
    } else {
        if (tx == 0) {
            if (node_data_out[i].left == best_i) {
                node_data_out[i].left = make_node(best_i, n, best*0.5, node_data_out, next_node_data_idx)->id;
            } else {
                node_data_out[i].right = make_node(best_i, n, best*0.5, node_data_out, next_node_data_idx)->id;
            }
        }
    }
    if (tx == 0) {
        node_data_out[best_i].dist = 0.5*best;
        nodes[n] = NULL;
    }
}

__global__ void parse_data(const char* D_raw, const uint32_t* D_raw_offsets, uint64_t n, uint32_t D_numlines, float* D_dst, uint8_t write_to_start) { // parse string to float on GPU - process each line separately
    uint32_t tx = threadIdx.x + blockDim.x*blockIdx.x;
    uint32_t stride = blockDim.x * gridDim.x;
    uint32_t i, j;
    char c;
    float f, m;
    for (i=tx; i<D_numlines; i+=stride) {
        float* D_dst_row = D_dst + (n+i)*(n+i-1)/2 - (write_to_start ? n*(n-1)/2 : 0);
        const char* D_raw_row = D_raw + D_raw_offsets[i] - 1;
        for (j=0; j<(n+i); j++) {
            f = 0;
            do {
                c = *(++D_raw_row);
            } while (c == ' ' || c == '\t');
            f = 0.0;
            while (c >= '0' && c <= '9') {
                f = 10.0*f + (float)(c - '0');
                c = *(++D_raw_row);
            }
            if (c == '.') {
                m = 1.0;
                c = *(++D_raw_row);
                while (c >= '0' && c <= '9') {
                    m *= 0.1;
                    f += (float)(c - '0') * m;
                    c = *(++D_raw_row);
                }
            }
            *(D_dst_row++) = f;
        }
    }
}

__global__ void nj_finalize(uint32_t N, float* D, Node** nodes, uint32_t* result_out) {
    uint32_t remnants[3];
    uint32_t i, j, k;
    double dist_i, dist_j, dist_k;
    j=0;
    for (i=0; i<N; i++) {
        if (nodes[i] != NULL) {
            remnants[j++] = i;
        }
    }
    i = remnants[0], j = remnants[1], k = remnants[2];
    dist_i = (D[index(j,i)] + D[index(k,i)] - D[index(k,j)]) * 0.5;
    dist_j = D[index(j,i)] - dist_i;
    dist_k = D[index(k,i)] - dist_i;
    nodes[i]->dist = dist_i;
    nodes[j]->dist = dist_j;
    nodes[k]->dist = dist_k;

    result_out[0] = nodes[i]->id;
    result_out[1] = nodes[j]->id;
    result_out[2] = nodes[k]->id;
}

/*  Compute r[i] =  Σ_k d(i,k) / (n-2).
    Then Q(i,j) = (n-2)*d(i,j) - Σ_k [d(i,k) + d(j,k)] becomes
    Q(i,j) = (n-2)*[d(i,j) - r[i] - r[j]], and because we only
    care about argmin[i,j] of Q, we can ignore the factor n-2.

    This can also be used for the distance to nodecalculation
    δ(f,u) = 1/2 * d(f,g) + Σ_k [d(f,k) - d(g,k)] / (2*(n-2)),
    which becomes 1/2 * [d(f,g) + r[f] - r[g]]. */
__global__ void compute_r(uint32_t N, const float* D, float* r, uint8_t normalize) {
    uint32_t i, j;
    uint32_t tx = threadIdx.x + blockDim.x*blockIdx.x;
    uint32_t stride = blockDim.x * gridDim.x;
    float ri;
    for (i=tx; i<N; i+=stride) {
        ri = 0.0;
        for (j=0; j<i; j++) ri += D[i*(i-1)/2 + j];
        for (j=i+1; j<N; j++) ri += D[j*(j-1)/2 + i];
        r[i] = normalize ? (ri / (N - 2.0)) : ri;
    }
}

float* d_D;
float* d_D_tmp; float* d_D_swap; // temporary array/pointer needed for multithreading compaction
float* d_r;
char* d_D_raw;
uint32_t* d_D_raw_offsets;
Node* d_node_data_out;
Node** d_nodes;
uint32_t* d_result_out;
uint32_t* d_next_node_data_idx;
min_data* d_best_ij;

void dispatch_parser(const char* D_raw, const uint32_t* D_raw_offsets, uint32_t n, uint32_t numlines, uint32_t data_size, uint8_t write_to_start) {
    cudaMemcpy(d_D_raw, D_raw, data_size*sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_D_raw_offsets, D_raw_offsets, numlines*sizeof(uint32_t), cudaMemcpyHostToDevice);
    parse_data<<<256, 256>>>(d_D_raw, d_D_raw_offsets, n, numlines, d_D, write_to_start);
}

void dispatch_setup(uint32_t N, uint32_t N_full) {
    cudaMalloc(&d_D, N*(N-1)/2*sizeof(float));
    cudaMalloc(&d_r, N_full*sizeof(float));
    cudaMalloc(&d_D_raw, (1<<28)*sizeof(char));
    cudaMalloc(&d_D_raw_offsets, N*sizeof(uint32_t));
    cudaMalloc(&d_node_data_out, 2*N_full*sizeof(Node));
    cudaMalloc(&d_nodes, N_full*sizeof(Node*));
    cudaMalloc(&d_result_out, 3*sizeof(uint32_t));
    cudaMalloc(&d_next_node_data_idx, sizeof(uint32_t));
    cudaMalloc(&d_best_ij, ARGMIN_NB*sizeof(min_data));

    nj_setup<<<1, 1>>>(N_full, d_node_data_out, d_nodes, d_next_node_data_idx);
}

void dispatch_skeleton(uint32_t N, uint32_t N_full) {
    compute_r<<<32, 256>>>(N, d_D, d_r, 1);
    cudaMalloc(&d_D_tmp, N*(N-1)/2*sizeof(float));
    uint32_t next_compaction = next_compaction_n(N);
    for (uint32_t n=N; n>3; n--) { // N is current matrix size (including dead elements), n is number of remaining nodes
        if (n == next_compaction ) {
            nj_compact<<<1, COMPACT_NT>>>(N, d_D, d_D_tmp, d_r, d_nodes);
            d_D_swap = d_D; d_D = d_D_tmp; d_D_tmp = d_D_swap; // swap pointers to new data location
            N = n;
            next_compaction = next_compaction_n(N);
        }
        nj_argmin<<<ARGMIN_NB, ARGMIN_NT>>>(N, n, d_D, d_r, d_best_ij);
        nj_update<<<1, UPDATE_NT>>>(N, n, d_D, d_r, d_node_data_out, d_nodes, d_next_node_data_idx, d_best_ij);
    }
}

void dispatch_add(const char* D_raw, const uint32_t* D_raw_offsets, uint32_t n, uint32_t numlines, uint32_t data_size) {
    dispatch_parser(D_raw, D_raw_offsets, n, numlines, data_size, 1);
    for (uint32_t i=0; i<numlines; i++) {
        nj_add<<<1, ADD_NT>>>(n+i, d_D + i*n + i*(i-1)/2, d_r, d_node_data_out, d_nodes, d_next_node_data_idx);
    }
}

uint32_t* dispatch_finalize(uint32_t N_full, Node** node_data_out) {
    nj_finalize<<<1, 1>>>(N_full, d_D, d_nodes, d_result_out);

    Node* node_data = (Node*)malloc(2*N_full*sizeof(Node)); // preallocate enough space for all nodes
    uint32_t* result = (uint32_t*)malloc(3*sizeof(uint32_t)); // indices of the three top-level nodes

    cudaMemcpy(node_data, d_node_data_out, 2*N_full*sizeof(Node), cudaMemcpyDeviceToHost);
    cudaMemcpy(result, d_result_out, 3*sizeof(uint32_t), cudaMemcpyDeviceToHost);

    cudaFree(d_D);
    cudaFree(d_D_tmp);
    cudaFree(d_r);
    cudaFree(d_node_data_out);
    cudaFree(d_nodes);
    cudaFree(d_result_out);
    cudaFree(d_next_node_data_idx);
    cudaFree(d_best_ij);
    cudaFree(d_D_raw);
    cudaFree(d_D_raw_offsets);

    cudaDeviceSynchronize();

    *node_data_out = node_data;
    return result;
}

