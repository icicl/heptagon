#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "hep.h"

double time_us() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ((double)ts.tv_sec * 1000000.0) + ((double)ts.tv_nsec / 1000.0);
}

int main(int argc, char** argv) {
    char* arg_in = NULL;
    uint8_t arg_flatten = 1;
    uint8_t arg_silent = 0;
    int32_t arg_skeleton_size = 20000;
    for (int argi=1; argi<argc; argi++) {
        if (strcmp(argv[argi], "--in") == 0) {
            if (++argi < argc) arg_in = argv[argi];
        } else if (strcmp(argv[argi], "--unflatten") == 0) {
            arg_flatten = 0;
        } else if (strcmp(argv[argi], "--silent") == 0) {
            arg_silent = 1;
        } else if (strcmp(argv[argi], "--skeleton-size") == 0) {
            if (++argi < argc) {
                arg_skeleton_size = atoi(argv[argi]);
                if (arg_skeleton_size <= 0) {
                    fprintf(stderr, "ERROR: Invalid skeleton size argument. Expected positive integer but got %s. Reverting to default value of 20000.\n", argv[argi]);
                    arg_skeleton_size = 20000;
                }
            }
        } else {
            fprintf(stderr, "Skipping unknown argument: %s\n", argv[argi]);
        }
    }
    if (arg_in == NULL) {
        fprintf(stderr, "Invalid invocation. Correct usage:\nheptagon --in filename [--flatten]\n");
        return 1;
    }

    uint32_t N, i, j, N_full;
    char* names;
    float timers[4] = {0};
    const int MAX_NAME_SIZE = 20;
    const int name_stride = MAX_NAME_SIZE + 1;
    char* line = NULL;
    uint64_t line_len = 0;

    timers[0] = -time_us();
    FILE* fp = fopen(arg_in, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error opening file %s.\n", arg_in);
        return 1;
    }
    if (getline(&line, &line_len, fp) == -1 || (N = N_full = atoi(line)) == 0) {
        fprintf(stderr, "ERROR: Expected phylip format distance matrix. The first line must be a decimal integer representing the number of taxa.\n");
        return 1;
    }
    fprintf(stderr, "Detected phylip matrix with N=%d.\n", N);
    if (N_full < 3) {
        fprintf(stderr, "Must supply matrix with at least 3 taxa. Exiting.\n");
        return 1;
    }
    if (N > arg_skeleton_size) {
        N = arg_skeleton_size;
        fprintf(stderr, "Provided matrix has %d taxa. Building skeleton from %d of them...\n", N_full, N);
    }
    dispatch_setup(N, N_full);
    names = malloc(N_full*name_stride*sizeof(char));
    const uint32_t D_raw_size = (1<<28);
    char* D_raw = malloc(D_raw_size*sizeof(char));
    uint32_t* D_raw_offsets = malloc(N*sizeof(uint32_t));
    int64_t len;
    uint32_t D_raw_offset = 0;
    uint32_t D_raw_cur_row = 0;

    fprintf(stderr, "Processing matrix file...\n");
    for (i=0; i<N; i++) {
        len = getline(&line, &line_len, fp);
        if (len == -1) {
            fprintf(stderr, "ERROR: Failed to read line #%d of %s.\n", i+2, arg_in);
            return 1;
        }
        if (D_raw_offset + len > D_raw_size) {
            dispatch_parser(D_raw, D_raw_offsets, i-D_raw_cur_row, D_raw_cur_row, D_raw_offset, 0);
            D_raw_offset = 0;
            D_raw_cur_row = 0;
        }

        for (j=0; line[j] != ' ' && line[j] != '\t'; j++) ;
        line[j] = '\0';
        if (j > MAX_NAME_SIZE) {
            fprintf(stderr, "WARNING: Max allowed taxon name length is %d. Truncating name %s.\n", MAX_NAME_SIZE, line);
            j = MAX_NAME_SIZE;
            line[j] = '\0';
        }
        memcpy(names + i*name_stride, line, j+1);
        memcpy(D_raw + D_raw_offset, line+j+1, len-j-1);
        D_raw_offsets[D_raw_cur_row] = D_raw_offset;
        D_raw_offset += (len-j-1);
        D_raw_cur_row++;
    }

    dispatch_parser(D_raw, D_raw_offsets, i-D_raw_cur_row, D_raw_cur_row, D_raw_offset, 0);
    timers[0] += time_us();
    fprintf(stderr, "Loaded data matrix %s in %.1lfs.\n", arg_in, (timers[0])/1000000.0);

    timers[1] -= time_us();
    dispatch_skeleton(N, N_full);
    timers[1] += time_us();
    fprintf(stderr, "Generated skeleton in %.1lfs.\n", (timers[1])/1000000.0);



    D_raw_offset = 0;
    D_raw_cur_row = 0;
    for (i=N; i<N_full; i++) {
        len = getline(&line, &line_len, fp);
        if (len == -1) {
            fprintf(stderr, "ERROR: Failed to read line #%d of %s.\n", i+2, arg_in);
            return 1;
        }
        if (D_raw_offset + len > D_raw_size || (i * D_raw_cur_row) > (N*(N-1)/2)) { // TODO simplify check
            dispatch_add(D_raw, D_raw_offsets, i-D_raw_cur_row, D_raw_cur_row, D_raw_offset);
            D_raw_offset = 0;
            D_raw_cur_row = 0;
        }

        for (j=0; line[j] != ' ' && line[j] != '\t'; j++) ;
        line[j] = '\0';
        if (j > MAX_NAME_SIZE) {
            fprintf(stderr, "WARNING: Max allowed taxon name length is %d. Truncating name %s.\n", MAX_NAME_SIZE, line);
            j = MAX_NAME_SIZE;
            line[j] = '\0';
        }
        memcpy(names + i*name_stride, line, j+1);
        memcpy(D_raw + D_raw_offset, line+j+1, len-j-1);
        D_raw_offsets[D_raw_cur_row] = D_raw_offset;
        D_raw_offset += (len-j-1);
        D_raw_cur_row++;
    }
    dispatch_add(D_raw, D_raw_offsets, i-D_raw_cur_row, D_raw_cur_row, D_raw_offset);


    fclose(fp);

    Node* node_data;
    uint32_t* results = dispatch_finalize(N_full, &node_data);
    if (!arg_silent) {
        if (arg_flatten) {
            print_node_flattened(node_data, results[0], names, name_stride, "(", ",");
            print_node_flattened(node_data, results[1], names, name_stride, "", ",");
            print_node_flattened(node_data, results[2], names, name_stride, "", ");\n");
        } else {
            print_node(node_data, results[0], names, name_stride, "(\n", ",\n");
            print_node(node_data, results[1], names, name_stride, "", ",\n");
            print_node(node_data, results[2], names, name_stride, "", ");\n");
        }
    }
    free(names);
    free(node_data);

}

