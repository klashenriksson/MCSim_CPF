#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>
#include <math.h>

#include "data.h"

int serialize_par(Par* par, char* buff)
{
  return sprintf(buff, "N %d nblock %d nsamp %d ntherm %d seed %d\n", par->N, par->nblock, par->nsamp, par->ntherm, par->seed);
}

void deserialize_par(FILE* f, Par* out_par)
{
  fscanf(f, "N %d nblock %d nsamp %d ntherm %d seed %d\n", &out_par->N, &out_par->nblock, &out_par->nsamp, &out_par->ntherm, &out_par->seed);
}

void datafile_get_filename(Par* par, char* filename_buff)
{
  #ifndef SURV_BIAS
  sprintf(filename_buff, "data/data_N%d.txt", par->N);
  #else
  sprintf(filename_buff, "data/data_surv_bias_N%d.txt", par->N);
  #endif
}

void datafile_write_block_results(FILE* datafile, result_t r, int block)
{
  fprintf(datafile, "block %d S2 %lf\n", block, r.S2);
}

int datafile_read_block_results(FILE* datafile, result_t* out_results)
{
    return fscanf(datafile, "block %*d S2 %lf\n", &out_results->S2);
}

int dir_exists(char* dirname)
{
  struct stat sb;
  return stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode);
}

void create_if_not_exists(char* dirname)
{
    if(!dir_exists(dirname))
    {
        #ifdef _WIN32
        _mkdir(dirname);
        #else
        mkdir(dirname, 0700);
        #endif
    }
}