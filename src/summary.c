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
  double variance_m;
  int nblocks;
} datafile_t;

datafile_t read_file(FILE* f, Par* par)
{
  datafile_t datafile;
  memset(&datafile, 0, sizeof(datafile));

  result_t r;
  double avg_m2 = 0.;
  
  int block = 0;
  while (datafile_read_block_results(f, &r, &block) != EOF) {
    datafile.avg_e += r.e;
    datafile.avg_c += r.c;
    datafile.avg_m += r.m;
    datafile.nblocks += 1;
    avg_m2 += r.m*r.m;
  }

  datafile.avg_e /= datafile.nblocks;
  datafile.avg_c /= datafile.nblocks;
  datafile.avg_m /= datafile.nblocks;
  avg_m2 /= datafile.nblocks;

  datafile.variance_m = 1.0/(datafile.nblocks-1) * (avg_m2 - datafile.avg_m*datafile.avg_m);

  return datafile;
}

// The result function should at least print out system size, temperature,
// magnetization, error bars for magnetization, and Binder's cumulant.
void result(Par *par, datafile_t* df)
{
  printf("L %3.3d T %8f m %8f m_ebars %8f\n", par->L, par->t, df->avg_m, df->variance_m);
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
  }
}
