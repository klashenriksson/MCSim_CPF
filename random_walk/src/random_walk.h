#ifndef RANDOM_WALK_H
#define RANDOM_WALK_H

#define FNAMESIZE 64

typedef struct result {
  double S2;
  double S2_Var;
} result_t;

typedef struct Par {
  int L;
  int N;
  int seed;
  int nblock;
  int nsamp;
} Par;

// In config.c
extern int write_config(Par *par, int *spin, char *fname);
extern int read_config(Par *par, int *spin, char *fname);

#endif