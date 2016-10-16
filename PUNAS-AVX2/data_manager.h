#ifndef DATA_MANAGER_H_
#define DATA_MANAGER_H_

#include<iostream>
#include<fstream>
#include<string.h>
#include<stdlib.h>
#include <ctype.h>
#include <omp.h>

#define D_ERROR_RATE 0.04
#define D_KMER_LEN 12
#define BIT_PARA 32
#define AVX_SIMD 8
#define MAX_REF 1000000000
#define MAX_READ 10000000
#define MAX_SINGLE_READ 1024
#define MAX_POS 200000000
#define SLICE_MAX_POS 10000000
typedef unsigned int WORD;

class data_manager{

public:

	/* input related */
	int ref_name_len;
    int ref_len;
    char* ref_name;
    char* ref;
	
	int read_num;
	int read_len;
	int read_name_len;

	char *read_content;
	char *read_name;
	
	/* Hash Table related */
    int window_size;
    float error_rate;
    int error_len;
	int pos_size ;
	int index_size;
	int* hv_pos;
	int* hv_index;

	/* bit  related */
	int ref_bit_len;
	WORD* refh;
	WORD* refl;
	
	int read_bit_len;
	WORD* readh;
	WORD* readl;

	/*position related*/
	int read_index;
	size_t total_candidates;
	int*  read_maps;
	int* pos_buffers;
	WORD* shd_scores;
	WORD* pos_hbits;
	WORD* pos_lbits;

	void init_data_manager(char* ref_file,char* read_file,float error_rate);
	
	void load_reads(char* read_file);
	void load_ref(char* ref_file);
	void generateHashTable();
	size_t hashVal(char *seq, int* flag);
};

#endif
