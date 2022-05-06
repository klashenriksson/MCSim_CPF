#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>

#include "data.h"

int serialize_par(Par* par, char* buff)
{
  return sprintf(buff, "L %3.3d T %5.3f nblock %3.3d nsamp %3.3d ntherm %3.3d seed %d\n", par->L, par->t, par->nblock, par->nsamp, par->ntherm, par->seed);
}

void deserialize_par(FILE* f, Par* out_par)
{
  fscanf(f, "L %d T %lf nblock %d nsamp %d ntherm %d seed %d\n", &out_par->L, &out_par->t, &out_par->nblock, &out_par->nsamp, &out_par->ntherm, &out_par->seed);
}

void datafile_get_filename(Par* par, char* filename_buff)
{
  sprintf(filename_buff, "data/data_L%3.3d_T%5.3f.txt", par->L, par->t);
}

void datafile_write_block_results(FILE* datafile, result_t r, int block)
{
  char buff[256] = {0};
  fprintf(datafile, "block %d e %8f c %8f m %8f m2 %8f m4 %8f\n", block, r.e, r.c, r.m, r.m2, r.m4);
}

int datafile_read_block_results(FILE* datafile, result_t* out_results, int* out_block)
{
    return fscanf(datafile, "block %d e %lf c %lf m %lf m2 %lf m4 %lf\n", out_block, &out_results->e, &out_results->c, &out_results->m, &out_results->m2, &out_results->m4);
}

void spincorrfile_get_filename(Par* par, char* filename_buff)
{
  sprintf(filename_buff, "corr/corr_L%d_T%5.3f.txt", par->L, par->t);
}

void spincorrfile_write_block_results(FILE* spincorrfile, double* correlations, int L, int block)
{
  for (int i = 0; i < L; i++)
  {
    fprintf(spincorrfile, "%8f\n", correlations[i]);
  }
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