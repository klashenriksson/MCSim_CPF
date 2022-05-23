#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include "percolation.h"
#include "data.h"

typedef struct datafile
{
  double avg_perc_prob;
  double avg_block_perc_prob2;
  int nblocks;
} datafile_t;

datafile_t read_file(FILE* f, Par* par)
{
  datafile_t datafile;
  memset(&datafile, 0, sizeof(datafile));

  result_t r;

  while (datafile_read_block_results(f, &r) != EOF) {
    datafile.avg_perc_prob += r.perc_prob;
    datafile.avg_block_perc_prob2 += r.perc_prob*r.perc_prob;
    datafile.nblocks += 1;
  }

  datafile.avg_perc_prob /= datafile.nblocks;
  datafile.avg_block_perc_prob2 /= datafile.nblocks;

  return datafile;
}

void result(Par *par, datafile_t* df)
{
  double perc_variance = 1.0/(fmax(df->nblocks-1,1)) * (df->avg_block_perc_prob2 - df->avg_perc_prob*df->avg_perc_prob);
  printf("L %d p %8f perc_prob %8f perc_prob_variance %8f\n", par->L, par->p, df->avg_perc_prob, perc_variance);
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
