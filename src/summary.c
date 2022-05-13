#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include "ising.h"
#include "data.h"

typedef struct datafile
{
  double avg_e;
  double avg_c;
  double avg_m;
  double avg_block_m2;
  double avg_m2;
  double avg_m4;
  double* correlations;
  int nblocks;
} datafile_t;

datafile_t read_file(FILE* f, Par* par)
{
  datafile_t datafile;
  memset(&datafile, 0, sizeof(datafile));

  result_t r;
  
  double* correlations = malloc(sizeof(double) * par->L);

  while (datafile_read_block(f, &r, correlations, par->L) != EOF) {
    datafile.avg_e += r.e;
    datafile.avg_c += r.c;
    datafile.avg_m += r.m;
    datafile.avg_block_m2 += r.m*r.m;
    datafile.avg_m2 += r.m2;
    datafile.avg_m4 += r.m4;
    datafile.nblocks += 1;
  }

  datafile.avg_e /= datafile.nblocks;
  datafile.avg_c /= datafile.nblocks;
  datafile.avg_m /= datafile.nblocks;
  datafile.avg_block_m2 /= datafile.nblocks;
  datafile.avg_m2 /= datafile.nblocks;
  datafile.avg_m4 /= datafile.nblocks;
  datafile.correlations = correlations;

  return datafile;
}

// The result function should at least print out system size, temperature,
// magnetization, error bars for magnetization, and Binder's cumulant.
void result(Par *par, datafile_t* df)
{
  double m_variance = 1.0/(fmax(df->nblocks-1,1)) * (df->avg_block_m2 - df->avg_m*df->avg_m);
  double Q = df->avg_m2*df->avg_m2/df->avg_m4;
  printf("L %d T %8f m %8f m_variance %8f Q %8f\n", par->L, par->t, df->avg_m, m_variance, Q);
}

void result_correlations(FILE* corrfile, Par *par, datafile_t* df)
{
    fprintf(corrfile, "L %d T %8f\nspin_corr:\n", par->L, par->t);
    for(int i = 0; i < par->L; i++)
    {
      fprintf(corrfile, "x %d value %8f\n", i, df->correlations[i]);
    }
}


int main(int argc, char *argv[])
{
  int i, nblock;
  char fname[FNAMESIZE];
  Par *par = malloc(sizeof(Par));

  while (fgets(fname, FNAMESIZE, stdin)) {
    // Change the trailing '\n' into a null character.
    char *s = strchr(fname, '\n');
    if (s - fname < FNAMESIZE)
      *s = '\0';

    // It is easier to use files with formatted data (written with fprintf).
    // FILE *stream;
    // stream = fopen(fname, "r");
    // if (!stream) {
    //   fprintf(stderr, "*** Can't open file %s.\n", fname);
    //   continue;
    // }

    FILE* f = fopen(fname, "r");
    if (!f)
    {
      fprintf(stderr, "*** Can't open file %s.\n", fname);
      continue;
    }
    deserialize_par(f, par);
    serialize_par(par, fname);
    datafile_t df = read_file(f, par);
    fclose(f);

    if (df.nblocks > 0)
      result(par, &df);

    //  Also print correlations
    char corrname[256] = {0};
    #ifndef TRI
    sprintf(corrname, "spincorr/tri_L%d_T%8f.txt", par->L, par->t);
    #else
    sprintf(corrname, "spincorr/tri_L%d_T%8f.txt", par->L, par->t);
    #endif

    create_if_not_exists("spincorr/");
    FILE* f_spincorr = fopen(corrname, "w");
    if (!f_spincorr)
    {
      fprintf(stderr, "*** Can't open spin correlation file %s\n.", corrname);
      continue;
    }

    result_correlations(f_spincorr, par, &df);
    fclose(f_spincorr);

    free(df.correlations);
  }
}
