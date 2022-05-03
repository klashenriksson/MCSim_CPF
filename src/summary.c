#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include "ising.h"
#include "data.h"

#define NVAR 5

int read_file(Par *par, double *vsum, double *v2sum, FILE* f)
{
  int i, nbl = 0;
  double v[NVAR];

  for (i = 0; i < NVAR; i++)
    vsum[i] = v2sum[i] = 0.0;

  result_t r;
  int block = 0;
  while (datafile_read_block_results(f, &r, &block) != EOF) {
  }

  return nbl;
}

// The result function should at least print out system size, temperature,
// magnetization, error bars for magnetization, and Binder's cumulant.
void result(Par *par, double *v, double *v2, int nblock)
{}


int main(int argc, char *argv[])
{
  int i, nblock;
  char fname[FNAMESIZE];
  double vsum[NVAR], v2sum[NVAR];
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
    nblock = read_file(par, vsum, v2sum, f);
    fclose(f);

    if (nblock)
      result(par, vsum, v2sum, nblock);
  }
}
