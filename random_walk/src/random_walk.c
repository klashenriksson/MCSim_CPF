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

int is_visited(int* walk_buff, int L, ivec2d_t pos)
{
  return walk_buff[wrap(pos.x, 0, L-1) + wrap(pos.y, 0, L-1)*L] == 1;
}

void set_visited(int* walk_buff, int L, ivec2d_t pos)
{
  walk_buff[wrap(pos.x, 0, L-1) + wrap(pos.y,0,L-1)*L] = 1;
}

void measure(Par *par, int* walk_buff, int squared_distance, double w, double* out_s2, double* out_s4)
{
  *out_s2 = w*squared_distance;
  *out_s4 = w*squared_distance*squared_distance;
}


result_t result(Par *par, double s2_sum, double w_s2_sum, double s4_sum, double w_sum, double w2_sum)
{
  double nsamps = par->nsamp;
  double s2_mean = s2_sum / nsamps;
  double s4_mean = s4_sum / nsamps;

  double s2_weighted_mean = w_s2_sum / w_sum;

  double unweighted_s2_variance = nsamps/(nsamps - 1) * (s4_mean - s2_mean*s2_mean);
  double b = w_sum*w_sum/w2_sum;

  double weighted_s2_variance = unweighted_s2_variance/b;

  printf("  %8f %8f %d\n", s2_weighted_mean, weighted_s2_variance, par->nsamp);

  result_t r;
  r.S2 = s2_weighted_mean;
  r.S2_Var = weighted_s2_variance;

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

#ifdef SURV_BIAS
int get_valid_dirs(int* walk_buff, ivec2d_t* rel_dirs, int L, ivec2d_t pos, ivec2d_t abs_dir)
{
  int valid_dirs = 0;
  for (int i = 0; i < n_rel_dirs; i++)
  {
    ivec2d_t rdir = relative_dirs[i];
    ivec2d_t dir = dir_rel_to_abs(abs_dir, rdir);
    ivec2d_t p = ivec_add(pos, dir);
    if (!is_visited(walk_buff, L, p))
    {
      rel_dirs[valid_dirs] = rdir;
      valid_dirs++;
    }
  }

  return valid_dirs;
}

int update(Par* par, int* walk_buff, double* out_w)
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

  const ivec2d_t origin = {par->L/2,par->L/2};

  // generate initial direction where all dirs all valid!
  unsigned int abs_r = uran() % n_abs_dirs;
  ivec2d_t dir = abs_dirs[abs_r];

  ivec2d_t old_pos = origin;
  ivec2d_t new_pos = ivec_add(old_pos, dir);

  set_visited(walk_buff, par->L, old_pos);
  set_visited(walk_buff, par->L, new_pos);

  int N = par->N;
  int steps = 1;
  double w = 1.f;
  while (steps < N)
  {
    unsigned int r = uran();
    ivec2d_t valid_rel_dirs[3] = {0};
    int n_valid_rel_dirs = get_valid_dirs(walk_buff, valid_rel_dirs, par->L, new_pos, dir);
    if (n_valid_rel_dirs == 0)
    {
      return -1;
    }
    ivec2d_t rel_dir = valid_rel_dirs[r % n_valid_rel_dirs];
    dir = dir_rel_to_abs(dir, rel_dir);
    double w_i = (double)n_valid_rel_dirs/(double)n_rel_dirs;
    w *= w_i;

    old_pos = new_pos;
    new_pos = ivec_add(new_pos, dir);

    set_visited(walk_buff, par->L, new_pos);
    steps++;
  }

  *out_w = w;
  return ivec_mag2(ivec_sub(origin, new_pos));
}
#else
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

  const ivec2d_t origin = {par->L/2,par->L/2};

  // generate initial direction where all dirs all valid!
  unsigned int abs_r = uran() % n_abs_dirs;
  ivec2d_t dir = abs_dirs[abs_r];

  ivec2d_t old_pos = origin;
  ivec2d_t new_pos = ivec_add(old_pos, dir);

  set_visited(walk_buff, par->L, old_pos);
  set_visited(walk_buff, par->L, new_pos);

  int N = par->N;
  int steps = 1;
  while (steps < N)
  {
    unsigned int r = uran();
    ivec2d_t rel_dir = relative_dirs[r % n_rel_dirs];
    ivec2d_t heading = dir;
    dir = dir_rel_to_abs(dir, rel_dir);

    old_pos = new_pos;
    new_pos = ivec_add(new_pos, dir);

    if (is_visited(walk_buff, par->L, new_pos)) {
        return -1;
    }

    set_visited(walk_buff, par->L, new_pos);
    steps++;
  }

  return ivec_mag2(ivec_sub(origin, new_pos));
}
#endif

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
  int i, iblock, isamp, istep;

  char datafilename[256] = {0};
  datafile_get_filename(par, datafilename);
  FILE* data_file = fopen(datafilename, "a");
  if(!data_file)
  {
    fprintf(stderr, "Error! Expected datafile to be created. Aborting.\n");
    return -1;
  }

  // *** Write out information about the run: size, temperature,...
  #ifndef SURV_BIAS
  printf("2D Lattice Random Walk model.\n");
  #else
  printf("2D Lattice Random Walk model (survival bias).\n");
  #endif

  printf("\n====    N = %d    ====\n", par->N);
  printf("\nnblock   nsamp   seed\n");
  printf("   %5d   %5d   %d\n", par->nblock, par->nsamp, par->seed);


  printf("\n  S2      S2_Var        nsamps used\n");

  for (iblock = 0; iblock < par->nblock; iblock++) {
    double block_w_tot = 0.f;
    double block_w2_tot = 0.f;
    double block_s2 = 0.f;
    double block_s4 = 0.f;
    double block_weighted_s2 = 0.f;

    for (isamp = 0; isamp < par->nsamp; isamp++) {
      double sample_w = 1.f;
      #ifdef SURV_BIAS
      int squared_distance = update(par, walk_buff, &sample_w);
      #else
      int squared_distance = update(par, walk_buff);
      #endif
      
      if (squared_distance == -1) {
        //discard this sample
        isamp--;
        continue;
      }

      double weighted_sample_s2 = 0.f;
      double weighted_sample_s4 = 0.f;
      measure(par, walk_buff, squared_distance, sample_w, &weighted_sample_s2, &weighted_sample_s4);

      block_s2 += weighted_sample_s2 / sample_w;
      block_s4 += weighted_sample_s4 / sample_w;
      block_weighted_s2 += weighted_sample_s2;
      block_w_tot += sample_w;
      block_w2_tot += sample_w*sample_w;
    }

    write_config(par, walk_buff, fname);
    result_t r = result(par, block_s2, block_weighted_s2, block_s4, block_w_tot, block_w2_tot);
    datafile_write_block_results(data_file, r, iblock);
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
    printf("Optional arguments (with defaults) nblock=%d nsamp=%d seed=%d\n",
	   par->nblock, par->nsamp, par->seed);
    exit(EXIT_SUCCESS);
  }

  // Interpret the commands given in the argument list.
  for (iarg = 1; iarg < argc; iarg++)
    if (!read_args(par, argv[iarg]))
      exit(EXIT_FAILURE);

  free(par);
  return 0;
}
