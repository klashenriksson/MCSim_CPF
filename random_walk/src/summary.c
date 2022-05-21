#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include "random_walk.h"
#include "data.h"

typedef struct datafile
{
  double avg_s2;
  double avg_block_s4;
  int nblocks;
} datafile_t;

datafile_t read_file(FILE* f, Par* par)
{
  datafile_t datafile;
  memset(&datafile, 0, sizeof(datafile));

  result_t r;

  while (datafile_read_block_results(f, &r) != EOF) {
    datafile.avg_s2 += r.S2;
    datafile.avg_block_s4 += r.S2*r.S2;
    datafile.nblocks += 1;
  }

  datafile.avg_s2 /= datafile.nblocks;
  datafile.avg_block_s4 /= datafile.nblocks;

  return datafile;
}

void result(Par *par, datafile_t* df)
{
  double s2_variance = 1.0/(fmax(df->nblocks-1,1)) * (df->avg_block_s4 - df->avg_s2*df->avg_s2);
  printf("N %d s2 %8f s2_variance\n", par->N, df->avg_s2, s2_variance);
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
