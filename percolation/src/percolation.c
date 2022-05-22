#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#include "ran.h"
#include "percolation.h"

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

ivec2d_t ivec_sub(ivec2d_t a, ivec2d_t b)
{
  return ivec_new(a.x - b.x, a.y - b.y);
}

int ivec_equal(ivec2d_t a, ivec2d_t b)
{
  return a.x == b.x && a.y == b.y ? 1 : 0;
}

// Global variables
char fname[FNAMESIZE];	// Name for config files

const ivec2d_t relative_dirs[] = {
    {1, 0},// forward
    {0, 1}, // right
    {0, -1}, // left
  };
const int n_rel_dirs = sizeof(relative_dirs)/sizeof(ivec2d_t);

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

void measure(Par *par, int* lattice, int squared_distance, double w, double* out_s2)
{
  *out_s2 = w*squared_distance;
}

result_t result(Par *par, double s2, double w_tot, int final)
{
  double s2_mean = s2 / w_tot;

  if (final)
  {
    printf("  --------  --------  --------\n");
  }
  printf("  %8f\n", s2_mean);

  result_t r;
  r.S2 = s2_mean;

  return r;
}

ivec2d_t dir_rel_to_abs(ivec2d_t curr_dir, ivec2d_t rel_dir)
{
  ivec2d_t dir;
  if (curr_dir.x == 0 && curr_dir.y == 1) { //heading downm
    dir.x = rel_dir.y;
    dir.y = rel_dir.x;
  } else if (curr_dir.x == 0 && curr_dir.y == -1) { //heading up
    dir.x = -rel_dir.y;
    dir.y = -rel_dir.x;
  } else if (curr_dir.x == 1 && curr_dir.y == 0) { //heading right
    dir.x = rel_dir.x;
    dir.y = rel_dir.y;
  } else if (curr_dir.x == -1 && curr_dir.y == 0) { //heading left
    dir.x = -rel_dir.x;
    dir.y = -rel_dir.y;
  }
  return dir;
}

int update(Par* par, int* lattice)
{
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

int mc(Par *par, int *lattice)
{
  int i, iblock, isamp, istep, ntherm = par->ntherm;

  // *** Read in the configuration for the present parameters if already present.
  if (read_config(par, lattice, fname))
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
  printf("2D Lattice Percolation model.\n");

  printf("\n====    L = %d    ====\n", par->L);
  printf("\nntherm  nblock   nsamp   seed\n");
  printf(" %5d   %5d   %5d   %d\n", ntherm, par->nblock, par->nsamp, par->seed);


  printf("\n  S2\n");

  for (iblock = 0; iblock < par->nblock; iblock++) {

    for (isamp = 0; isamp < par->nsamp; isamp++) {
      double sample_s2 = 0.f;
      update(par, lattice);
      
      //measure(par, walk_buff, squared_distance, sample_w, &sample_s2);
    }

    write_config(par, lattice, fname);
    //result_t r = result(par, block_s2, tot_w, 0);
    //datafile_write_block_results(data_file, r, iblock);
  }

  fclose(data_file);

  return 1;
}

int initialize_mc(Par *par, int *lattice) {
  char *f2;

  if (!par->L) {
    printf("Give walk length N!\n");
    return 0;
  }


  init_ran(par->seed);

  sprintf(fname, "%d", par->L);

  check_data_file(par);

  return mc(par, lattice);
}


int read_args(Par *par, char *arg)
{
  int st;
  static int *lattice = NULL;
  char *s;

  if (!strcmp(arg, "run"))
    return initialize_mc(par, lattice);

  s = strchr(arg, '=');

  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'.", arg);
    return 0;
  }

  *s++ = '\0';

  if (!strcmp(arg, "N")) {
    par->L = strtod(s, NULL);
    lattice = realloc(lattice, par->L * par->L * sizeof(int));
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
    return read_config(par, lattice, s);
  }


  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}

int main(int argc, char *argv[])
{
  int i, iarg;
  Par *par = malloc(sizeof(Par));

  par->L = 0;
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
