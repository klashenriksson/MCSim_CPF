#ifndef DATA_H
#define DATA_H

#include "ising.h"

int serialize_par(Par* par, char* buff);
void deserialize_par(FILE* f, Par* out_par);

void datafile_get_filename(Par* par, char* filename_buff);
void datafile_write_block_results(FILE* datafile, result_t r, int block);
int datafile_read_block_results(FILE* datafile, result_t* out_results, int* out_block);
int dir_exists(char* dirname);
int create_if_not_exists(char* dirname);

#endif