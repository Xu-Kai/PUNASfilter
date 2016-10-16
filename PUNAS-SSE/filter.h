#ifndef NODE_COMPUTE_H_
#define NODE_COMPUTE_H_

#include <memory.h>
#include <stdint.h>
#include <omp.h>
#include <sys/sysinfo.h>
#include <sys/time.h>

#include <emmintrin.h>
#include <tmmintrin.h>

#include "data_manager.h"

int filter_reads(data_manager* dm,int SIMD_WIDTH);
int sse_shd(data_manager* dm,int SIMD_WIDTH,int pos_nums);
int find_candidates(data_manager* dm,int SIMD_WIDTH);
int find_read_pos(data_manager* dm, int read_id, int* kmer_index,int* kmer_value,int SIMD_WIDTH,int& pos_index);
void fill_regions(data_manager* dm, int pos_id,int SIMD_WIDTH);
void flip_false_zero(__m128i ham_mask, __m128i pre_ham, __m128i& fix_ham_mask);
void pop_count(__m128i& finalBV);
void flip_false_zero2(__m128i ham_mask, __m128i pre_ham, __m128i& fix_ham_mask);
void flip_false_zero3(__m128i ham_mask, __m128i pre_ham, __m128i& fix_ham_mask);


#endif
