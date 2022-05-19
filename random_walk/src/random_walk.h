#ifndef RANDOM_WALK_H
#define RANDOM_WALK_H

#define FNAMESIZE 64

typedef struct result {

} result_t;

typedef struct Par {
  int N;
  int seed;
  int ntherm;
  int nblock;
  int nsamp;
} Par;

// In config.c
extern int write_config(Par *par, int *spin, char *fname);
extern int read_config(Par *par, int *spin, char *fname);

#endif