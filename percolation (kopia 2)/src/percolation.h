#ifndef PERCOLATION_H
#define PERCOLATION_H

#define FNAMESIZE 64
#define D 2

typedef struct result {
  double perc_prob;
} result_t;

typedef struct Par {
  int L;
  double p;
  int seed;
  int ntherm;
  int nblock;
  int nsamp;
} Par;

#endif