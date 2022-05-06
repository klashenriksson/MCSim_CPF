#ifndef DATA_H
#define DATA_H

#include "ising.h"

int serialize_par(Par* par, char* buff);
void deserialize_par(FILE* f, Par* out_par);

void datafile_get_filename(Par* par, char* filename_buff);
void datafile_write_block_results(FILE* datafile, result_t r, int block);
int datafile_read_block_results(FILE* datafile, result_t* out_results, int* out_block);

void spincorrfile_get_filename(Par* par, char* filename_buff);
void spincorrfile_write_block_results(FILE* spincorrfile, double* correlations, int L, int block);
int spincorrfile_read_block_results(FILE* spincorrfile, double** out_correlations, int* outL);

int dir_exists(char* dirname);
void create_if_not_exists(char* dirname);

#endif