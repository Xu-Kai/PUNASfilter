#include "data_manager.h"
#include "tbb/parallel_sort.h"
#ifdef __GNUC__
#include <mm_malloc.h>
#endif

typedef struct {
	int key;
	int value;
} hv_loc;

void data_manager::init_data_manager(char* ref_file,char* read_file,float error_rate){

	load_reads(read_file);
	load_ref(ref_file);

	// specefic  args
	error_len = (int) ((float) read_len * error_rate);
	window_size = D_KMER_LEN;
	printf("error len %d, kmer len %d\n",error_len,window_size);
	
	generateHashTable();
	
	/*position related*/
	read_index=0;
	total_candidates=0;
	read_maps=(int*)malloc(MAX_POS*sizeof(int));
	pos_buffers=(int*)malloc(MAX_POS*sizeof(int));
	shd_scores=(WORD*)_mm_malloc(MAX_POS*sizeof(int),32);
	pos_hbits=(WORD*)_mm_malloc(MAX_POS*read_bit_len*sizeof(WORD),64);
	pos_lbits=(WORD*)_mm_malloc(MAX_POS*read_bit_len*sizeof(WORD),64);
	
}

void data_manager::load_ref(char* ref_file){
	
	printf("load ref\n");
	FILE *fp;
	char* line=NULL;
	size_t len=0;
	ssize_t read;
	
	ref_name_len=100;
	ref_len=0;
	ref_name=(char*)malloc(ref_name_len*sizeof(char));
	ref=(char*)malloc(MAX_REF*sizeof(char));
	
	fp=fopen(ref_file,"r");
	read=getline(&line,&len,fp);
	memcpy(ref_name,line,read-1);

	while(read!=-1){
		
		read=getline(&line,&len,fp);
		
		if(read==-1)
			break;
		
		if(read-1>0){
			memcpy(ref+ref_len,line,read-1);
			ref_len+=read-1;
		}
		
	}

	
	ref_bit_len=(ref_len+BIT_PARA-1)/BIT_PARA;
	refh=(WORD*)malloc(ref_bit_len*sizeof(WORD));
	refl=(WORD*)malloc(ref_bit_len*sizeof(WORD));	

	
	#pragma omp parallel for
	for (int i = 0; i < ref_len; i += BIT_PARA) {

		for (int j = 0; j < BIT_PARA; j++) {
			char c = toupper(ref[i + j]);
			switch (c) {
			case 'A':
				refh[i / BIT_PARA] = (refh[i / BIT_PARA] << 1) + 0;
				refl[i / BIT_PARA] = (refl[i / BIT_PARA] << 1) + 0;
				break;
			case 'C':
				refh[i / BIT_PARA] = (refh[i / BIT_PARA] << 1) + 1;
				refl[i / BIT_PARA] = (refl[i / BIT_PARA] << 1) + 1;
				break;
			case 'G':
				refh[i / BIT_PARA] = (refh[i / BIT_PARA] << 1) + 0;
				refl[i / BIT_PARA] = (refl[i / BIT_PARA] << 1) + 1;
				break;
			case 'T':
				refh[i / BIT_PARA] = (refh[i / BIT_PARA] << 1) + 1;
				refl[i / BIT_PARA] = (refl[i / BIT_PARA] << 1) + 0;
				break;
			default:
				refh[i / BIT_PARA] = (refh[i / BIT_PARA] << 1) + 0;
				refl[i / BIT_PARA] = (refl[i / BIT_PARA] << 1) + 0;
			}
		}

	}
	
	printf("REF--> %s with %d length %d bit length\n",ref_name,ref_len,ref_bit_len);
}

void data_manager::load_reads(char* read_file){
	printf("load reads\n");
	FILE *fp;
	char* line=NULL;
	size_t len=0;
	ssize_t read;
	
	fp=fopen(read_file,"r");
	read=getline(&line,&len,fp);
	read=getline(&line,&len,fp);
	read_len=read-1;
	
	
	read_name_len=100;
	read_num=0;
	read_name=(char*)malloc(MAX_READ*read_name_len*sizeof(char));
	long long int read_content_size=(long long int)(MAX_READ)*(long long int)(read_len);
	read_content=(char*)malloc(read_content_size*sizeof(char));
	read_bit_len=(read_len+BIT_PARA-1)/BIT_PARA;
	readh=(WORD*)malloc(MAX_READ*read_bit_len*sizeof(WORD));
	readl=(WORD*)malloc(MAX_READ*read_bit_len*sizeof(WORD));
	
	fseek(fp, 0, SEEK_SET);
	read=getline(&line,&len,fp);
	int count=0;
	while(read!=-1){
		read=getline(&line,&len,fp);
		count++;
	
		if(count%2==1){
			long long int offset=(long long int)(read_num)*(long long int)(read_len);
			memcpy(read_content+offset,line,read_len);
			read_num++;
		}
	}
	
	//找到read对应翻转之后的

#pragma omp parallel for
	for (int i = 0; i < read_num; i++) {
		for (int j = 0; j < read_len; j++) {
			char c = read_content[i * read_len + read_len - j - 1];
			switch (c) {
			case 'A':
				read_content[(read_num + i) * read_len + j] = 'T';
				break;
			case 'C':
				read_content[(read_num + i) * read_len + j] = 'G';
				break;
			case 'G':
				read_content[(read_num + i) * read_len + j] = 'C';
				break;
			case 'T':
				read_content[(read_num + i) * read_len + j] = 'A';
				break;
			default:
				read_content[(read_num + i) * read_len + j] = c;
			}
		}
	}
	read_num += read_num;

	
	
	int fixed_read_len = read_bit_len * BIT_PARA;
#pragma omp parallel for
	for (int i = 0; i < read_num; i++) {
		for (int j = 0; j < fixed_read_len; j++) {

			int add = i * (fixed_read_len / BIT_PARA) + j / BIT_PARA;

			if (j < read_len) {
				char c = read_content[i * read_len + j];
				switch (c) {
				case 'A':
					readh[add] = (readh[add] << 1) + 0;
					readl[add] = (readl[add] << 1) + 0;
					break;
				case 'C':
					readh[add] = (readh[add] << 1) + 1;
					readl[add] = (readl[add] << 1) + 1;
					break;
				case 'G':
					readh[add] = (readh[add] << 1) + 0;
					readl[add] = (readl[add] << 1) + 1;
					break;
				case 'T':
					readh[add] = (readh[add] << 1) + 1;
					readl[add] = (readl[add] << 1) + 0;
					break;
				default: // symbol==N
					readh[add] = (readh[add] << 1) + 0;
					readl[add] = (readl[add] << 1) + 0;
					break;
				}
			} else { //padding with A
				readh[add] = (readh[add] << 1) + 0;
				readl[add] = (readl[add] << 1) + 0;
			}

		}
	}
	
	
	printf("totally %d read with %d average length and %d bit length\n",read_num,read_len,read_bit_len);
	
	
}


bool compare(const hv_loc &a, const hv_loc &b) {
	if (a.key < b.key)
		return true;
	else if (a.key == b.key) {
		if (a.value < b.value)
			return true;
		else
			return false;
	} else
		return false;
}

size_t data_manager::hashVal(char *seq, int* flag) {
	int i = 0;
	size_t val = 0, numericVal = 0;
	*flag = 1;
	while (i < window_size) {
		switch (seq[i]) {
		case 'A':
			numericVal = 0;
			break;
		case 'C':
			numericVal = 1;
			break;
		case 'G':
			numericVal = 2;
			break;
		case 'T':
			numericVal = 3;
			break;
		default: {
			*flag = 0;
			return 0;
			break;
		}
		}
		val = (val << 2) | numericVal;
		i++;
	}
	return val;
}


void data_manager::generateHashTable() {
	printf("generate hash index \n");
	
	pos_size = ref_len - window_size + 1;
	index_size = (1 << (2 * window_size)) + 1;

	hv_loc* hv_loc_arr = (hv_loc*) malloc(pos_size * sizeof(hv_loc));
	hv_pos = (int*) malloc(pos_size * sizeof(int));
	hv_index = (int*) malloc(index_size * sizeof(int));
	memset(hv_index, 0, index_size * sizeof(int));

	int MAX_KEY = (1 << (2 * window_size)) + 10;
#pragma omp parallel for
	for (int i = 0; i < pos_size; i++) {
		size_t hv;
		int flag;
		hv = hashVal(ref + i, &flag);
		if (flag) {
			hv_loc_arr[i].key = hv;
			hv_loc_arr[i].value = i;
		} else {
			hv_loc_arr[i].key = MAX_KEY;
			hv_loc_arr[i].value = i;
		}
	}

	tbb::parallel_sort(hv_loc_arr, hv_loc_arr + pos_size, compare);

	int wrong = 0;

	for (int i = 0; i < pos_size; i++) {
		int cur_key = hv_loc_arr[i].key;
		if (cur_key < MAX_KEY) {
			hv_index[cur_key + 1]++;
		} else {
			wrong++;
		}
	}

	int max_kmer = 0;
	int mpos = 0;
	for (int i = 0; i < pos_size; i++) {
		int cur_key = hv_loc_arr[i].key;
		if (cur_key < MAX_KEY) {
			if (max_kmer < hv_index[cur_key + 1]) {
				max_kmer = hv_index[cur_key + 1];
				mpos = i;
			}
		}
	}

	for (int i = 0; i < index_size - 1; i++) {
		hv_index[i + 1] = hv_index[i] + hv_index[i + 1];
	}

#pragma omp parallel for
	for (int i = 0; i < pos_size; i++) {
		hv_pos[i] = hv_loc_arr[i].value;
	}

	free(hv_loc_arr);

}


