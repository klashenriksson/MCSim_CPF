#ifndef ISING_H
#define ISING_H

#include "int_queue.h"

#ifdef SSFAST
#define SS
#endif

#define FNAMESIZE 64

typedef struct result {
  double e;
  double m;
  double c;
} result_t;

typedef struct Par {
  int L;
  double t;
  int ntherm, nblock, nsamp, seed;
  #ifdef CLU
  int_queue_t queue;
  #endif
} Par;

// In config.c
extern int write_config(Par *par, int *spin, char *fname);
extern int read_config(Par *par, int *spin, char *fname);

#endif