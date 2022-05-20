#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#include "ran.h"
#include "random_walk.h"

#include "data.h"

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

typedef struct ivec2d {
  int x;
  int y;
} ivec2d_t;

int ivec_dot(ivec2d_t a, ivec2d_t b)
{
  return a.x * b.x + a.y * b.y;
}

int ivec_mag2(ivec2d_t a)
{
  return ivec_dot(a,a);
}

ivec2d_t ivec_new(int x, int y)
{
  ivec2d_t r;
  r.x = x;
  r.y = y;
  return r;
}

ivec2d_t ivec_add(ivec2d_t a, ivec2d_t b)
{
  return ivec_new(a.x + b.x, a.y + b.y);
}

int ivec_equal(ivec2d_t a, ivec2d_t b)
{
  return a.x == b.x && a.y == b.y ? 1 : 0;
}

// Global variables
char fname[FNAMESIZE];	// Name for config files

int wrap(int x, int min, int max)
{
  if (x < min)
  {
    int diff = x-min;
    return wrap(max + diff + 1, min, max);
  } else if (x > max)
  {
    int diff = x - max;
    return wrap(min + diff - 1, min, max);
  }

  return x;
}

void measure(Par *par, int* walk_buff, double* out_s2)
{
  *out_s2 = 0.f;

  int D_x = 0;
  int D_y = 0;
  for (int i = 0; i < par->N; i++)
  {
    int x_i = walk_buff[i * 2];
    int y_i = walk_buff[i * 2 + 1];
    D_x += x_i;
    D_y += y_i;
  }

  *out_s2 = D_x*D_x + D_y*D_y;
}


result_t result(Par *par, double s2, int ntot, int final)
{
  double s2_mean = s2 / ntot;

  if (final)
  {
    printf("  --------  --------  --------\n");
  }
  printf("  %8f\n", s2_mean);

  result_t r;
  r.S2 = s2_mean;

  return r;
}

int is_visited(int* walk_buff, int N, ivec2d_t pos)
{
  ivec2d_t p = ivec_new(0,0);
  for (int i = 0; i < N; i++)
  {

    if (ivec_equal(p, pos)) {
      return 1;
    }

    int d_x = walk_buff[i*2];
    int d_y = walk_buff[i * 2 + 1];

    p = ivec_add(p, ivec_new(d_x, d_y));
  }

  return 0;
}

int update(Par* par, int* walk_buff)
{
  const ivec2d_t abs_dirs[] = {
    {1, 0},
    {-1, 0},
    {0, 1},
    {0, -1}
  };
  const int n_abs_dirs = sizeof(abs_dirs)/sizeof(ivec2d_t);

  const ivec2d_t relative_dirs[] = {
    {1, 0},// forward
    {0, 1}, // right
    {0, -1}, // left
  };
  const int n_rel_dirs = sizeof(relative_dirs)/sizeof(ivec2d_t);

  // generate initial direction where all dirs all valid!
  unsigned int abs_r = uran() % n_abs_dirs;
  ivec2d_t dir = abs_dirs[abs_r];

  walk_buff[0] = dir.x;
  walk_buff[1] = dir.y;

  ivec2d_t old_pos = ivec_new(0,0);
  ivec2d_t new_pos = dir;

  for (int i = 1; i < par->N; i++)
  {
    unsigned int rel_idx = uran() % n_rel_dirs;
    ivec2d_t rel_dir = relative_dirs[rel_idx];
    ivec2d_t heading = dir;

    if (dir.x == 0 && dir.y == 1) { //heading up
      dir.x = rel_dir.y;
      dir.y = rel_dir.x;
    } else if (dir.x == 0 && dir.y == -1) { //heading down
      dir.x = -rel_dir.y;
      dir.y = -rel_dir.x;
    } else if (dir.x == 1 && dir.y == 0) { //heading right
      dir.x = rel_dir.x;
      dir.y = -rel_dir.y;
    } else if (dir.x== -1 && dir.y == 0) { //heading left
      dir.x = -rel_dir.x;
      dir.y = rel_dir.y;
    }

    old_pos = new_pos;
    new_pos = ivec_add(new_pos, dir);
    if (is_visited(walk_buff, i-1, new_pos)) {
      return 0;
    }

    walk_buff[i*2] = dir.x;
    walk_buff[i*2 + 1] = dir.y;
  }

  return 0;
}

void fwrite_par(FILE* f, Par* par) 
{
  char par_buff[256] = {0};
  int par_length = serialize_par(par, par_buff);
  fwrite(par_buff, sizeof(char), par_length, f);
}

void check_data_file(Par* par) {
  create_if_not_exists("data/");

  char filename[256] = {0};
  datafile_get_filename(par, filename);
  FILE* f = fopen(filename, "r");
  if (!f)
  {
    f = fopen(filename, "w");
    fwrite_par(f, par);
    fclose(f);
  }
}

int mc(Par *par, int *walk_buff)
{
  int i, iblock, isamp, istep, ntherm = par->ntherm;
  double acc, accept = 0.0;


  // *** Read in the configuration for the present parameters if already present.
  if (read_config(par, walk_buff, fname))
    ntherm = 0;

  char datafilename[256] = {0};
  datafile_get_filename(par, datafilename);
  FILE* data_file = fopen(datafilename, "a");
  if(!data_file)
  {
    fprintf(stderr, "Error! Expected datafile to be created. Aborting.\n");
    return -1;
  }

  // *** Write out information about the run: size, temperature,...
  printf("2D Lattice Random Walk modelm.\n");

  printf("\n====    N = %d    ====\n", par->N);
  printf("\nntherm  nblock   nsamp   seed\n");
  printf(" %5d   %5d   %5d   %d\n", ntherm, par->nblock, par->nsamp, par->seed);

  
  printf("\n  S2\n");

  #ifdef CORR
  corr_t corr = corr_create(par->nblock * par->nsamp);
  #endif

  // Thermalize the system 
  for (i = 0; i < ntherm; i++)
    update(par, walk_buff);

  double s2 = 0.f;
  for (iblock = 0; iblock < par->nblock; iblock++) {
    double block_s2 = 0.f;
    for (isamp = 0; isamp < par->nsamp; isamp++) {
      double sample_s2 = 0.f;
      accept += update(par, walk_buff);
      measure(par, walk_buff, &sample_s2);

      block_s2 += sample_s2;
    }

    write_config(par, walk_buff, fname);
    result_t r = result(par, block_s2, par->nsamp, 0);
    datafile_write_block_results(data_file, r, iblock);

    s2 += block_s2;
  }
  result(par, s2, par->nblock * par->nsamp, 1);

  fclose(data_file);

  acc = 0;
  printf("\nAcceptance: %5.2f\n", acc);

  return 1;
}

int initialize_mc(Par *par, int *walk_buff) {
  char *f2;

  if (!par->N) {
    printf("Give walk length N!\n");
    return 0;
  }


  init_ran(par->seed);

  sprintf(fname, "%d", par->N);

  check_data_file(par);

  return mc(par, walk_buff);
}


int read_args(Par *par, char *arg)
{
  int st;
  static int *walk_buff = NULL;
  char *s;

  if (!strcmp(arg, "run"))
    return initialize_mc(par, walk_buff);

  s = strchr(arg, '=');

  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'.", arg);
    return 0;
  }

  *s++ = '\0';

  if (!strcmp(arg, "N")) {
    par->N = strtod(s, NULL);
    walk_buff = realloc(walk_buff, par->N * sizeof(int) * 2);
    par->ntherm = 0;

    return 1;
  }

  if (!strcmp(arg, "ntherm")) {
    par->ntherm = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "nblock")) {
    par->nblock = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "nsamp")) {
    par->nsamp = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "seed")) {
    par->seed = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "read")) {
    return read_config(par, walk_buff, s);
  }


  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}

int main(int argc, char *argv[])
{
  int i, iarg;
  Par *par = malloc(sizeof(Par));

  par->N = 0;
  par->nblock = 1;
  par->nsamp = 10000;

  if (argc == 1) {
    printf("Usage: %s N=16\n", argv[0]);
    printf("Optional arguments (with defaults) nblock=%d nsamp=%d ntherm=%d seed=%d\n",
	   par->nblock, par->nsamp, par->ntherm, par->seed);
    exit(EXIT_SUCCESS);
  }

  // Interpret the commands given in the argument list.
  for (iarg = 1; iarg < argc; iarg++)
    if (!read_args(par, argv[iarg]))
      exit(EXIT_FAILURE);

  free(par);
  return 0;
}
