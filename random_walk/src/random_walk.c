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

int is_visited(int* walk_buff, int L, ivec2d_t pos)
{
  return walk_buff[wrap(pos.x, 0, L-1) + wrap(pos.y, 0, L-1)*L] == 1;
}

void set_visited(int* walk_buff, int L, ivec2d_t pos)
{
  walk_buff[wrap(pos.x, 0, L-1) + wrap(pos.y,0,L-1)*L] = 1;
}

void measure(Par *par, int* walk_buff, int steps, double* out_s2)
{
  *out_s2 = steps;
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

int update(Par* par, int* walk_buff)
{
  int L2 = par->L * par->L;
  memset(walk_buff, 0, L2 * sizeof(int));
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

  ivec2d_t old_pos = ivec_new(par->L/2,par->L/2);
  ivec2d_t new_pos = ivec_add(old_pos, dir);

  set_visited(walk_buff, par->L, old_pos);
  set_visited(walk_buff, par->L, new_pos);

  int N = par->N;
  int steps = 1;
  while (steps < N)
  {
    unsigned int rel_idx = uran() % n_rel_dirs;
    ivec2d_t rel_dir = relative_dirs[rel_idx];
    ivec2d_t heading = dir;
    dir = dir_rel_to_abs(dir, rel_dir);

    old_pos = new_pos;
    new_pos = ivec_add(new_pos, dir);
    if (is_visited(walk_buff, par->L, new_pos)) {
        return -1;
    }

    set_visited(walk_buff, par->L, new_pos);
    steps++;

    ivec2d_t nbor_0 = ivec_add(new_pos, abs_dirs[0]);
    ivec2d_t nbor_1 = ivec_add(new_pos, abs_dirs[1]);
    ivec2d_t nbor_2 = ivec_add(new_pos, abs_dirs[2]);
    ivec2d_t nbor_3 = ivec_add(new_pos, abs_dirs[3]);
    if (is_visited(walk_buff, par->L, nbor_0) && 
        is_visited(walk_buff, par->L, nbor_1) &&
        is_visited(walk_buff, par->L, nbor_2) &&
        is_visited(walk_buff, par->L, nbor_3))
    {
      //hard stuck
      /*for (int i = 0; i < par->L*par->L; i++)
      {
        int x = i % par->L;
        int y = i / par->L;
        ivec2d_t p = ivec_new(x,y);
        if (is_visited(walk_buff, par->L, p))
        {
          printf("v: %d, %d\n", p.x, p.y);
        }
      }
      printf("-- %d\n", steps);*/
      return steps;
    }
  }

  return steps;
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
  printf("2D Lattice Random Walk model.\n");

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
    int samples = par->nsamp;
    double block_s2 = 0.f;
    for (isamp = 0; isamp < par->nsamp; isamp++) {
      double sample_s2 = 0.f;
      int steps = update(par, walk_buff);
      if (steps == -1) {
        //discard this sample
        samples--;
        continue;
      }
      measure(par, walk_buff, steps, &sample_s2);

      block_s2 += sample_s2;
    }

    printf("nsamps used: %d\n", samples);
    write_config(par, walk_buff, fname);
    result_t r = result(par, block_s2, samples, 0);
    datafile_write_block_results(data_file, r, iblock);

    s2 += block_s2;
  }

  fclose(data_file);

  return 1;
}

int initialize_mc(Par *par, int *walk_buff) {
  char *f2;

  if (!par->L) {
    printf("Give walk length N!\n");
    return 0;
  }


  init_ran(par->seed);

  sprintf(fname, "%d", par->L);

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
    par->L = par->N*2 + 1;
    walk_buff = realloc(walk_buff, par->L * par->L * sizeof(int));
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

  par->L = 0;
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
