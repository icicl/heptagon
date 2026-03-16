#include <stdint.h>
#include "node.h"

#ifdef __cplusplus
extern "C" {
#endif
uint32_t compact(uint32_t N, float** D, Node** nodes, float* r);
void dispatch_setup(uint32_t N, uint32_t N_full);
void dispatch_parser(const char* D_raw, const uint32_t* D_raw_offsets, uint32_t n, uint32_t numlines, uint32_t data_size, uint8_t write_to_start);
void dispatch_skeleton(uint32_t N, uint32_t N_full);
void dispatch_add(const char* D_raw, const uint32_t* D_raw_offsets, uint32_t n, uint32_t numlines, uint32_t data_size);
uint32_t* dispatch_finalize(uint32_t N_full, Node** node_data_out);
#ifdef __cplusplus
}
#endif
