#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

#include "data_manager.h"
#include "filter.h"

int main(int argc, char* argv[]) {

	char* ref_file;
	char* read_file;
	float error_rate;

		if (argc >3) {
			ref_file = argv[1];
			read_file = argv[2];
			error_rate=atof(argv[3]);
		}else{
			printf(" Program needs at least two arugments: reference file name , read file name, error rate \n");
			printf(" For example :  ./eshd GRch37.fasta read.fasta error_rate\n");
			return 0;
		}

		printf("\n\n");
		printf(
				"*********************ESHD: 2016 shandong university hpc lab*****************\n");
		printf(
				"*******ESHD: Accelerating seed verification for next-generation sequencing read alignmenton multi-core and many-core architectures**********\n");
		printf("*******Error rate %f  kmer length %d **********\n",
		error_rate, D_KMER_LEN);
		printf("\n\n");


	data_manager dm;
	dm.init_data_manager(ref_file, read_file,error_rate);
	
	filter_reads(&dm,SSE_SIMD);
	
	return 0;

}

