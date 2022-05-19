#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#include "ran.h"
#include "random_walk.h"


#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif


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

void measure(Par *par)
{

}


result_t result(Par *par, int ntot, int final)
{
  if (final)
  {
    printf("  --------  --------  --------\n");
  }

  result_t r;

  return r;
}

int update(Par* par, int* spin, int draw)
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

int mc(Par *par, int *spin)
{
  int i, iblock, isamp, istep, ntherm = par->ntherm;
  double t = par->t, acc, accept = 0.0, L2 = par->L * par->L;
  double energy = 0.0, energy_sqrd = 0.0, magnetization = 0.0, mag2 = 0.0, mag4 = 0.0;


  // *** Read in the configuration for the present parameters if already present.
  if (read_config(par, spin, fname))
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
#ifdef CLU
  printf("2D Ising model with the Wolff cluster algorithm.\n");
  par->queue = int_queue_create(L2);
#else
  printf("2D Ising model with the Metropolis algorithm.\n");
#endif

#ifdef TRI
  printf("Using Triangular lattice.\n");
#endif

  printf("\n====    %d x %d     T = %g    ====\n", par->L, par->L, par->t);
  printf("\nntherm  nblock   nsamp   seed\n");
  printf(" %5d   %5d   %5d   %d\n", ntherm, par->nblock, par->nsamp, par->seed);

  
  printf("\n energy      cv        magn     \n");

  #ifdef CORR
  corr_t corr = corr_create(par->nblock * par->nsamp);
  #endif

  // Thermalize the system 
  for (i = 0; i < ntherm; i++)
    update(par, spin, 0);

  #ifdef DRAW_GRAPHICS
  int num_imgs = 5;
  int total_iterations = par->nblock * par->nsamp;
  int iteration_quotient = total_iterations/num_imgs; 
  #endif

  double* block_spin_corrs = malloc(sizeof(double)*par->L);

  for (iblock = 0; iblock < par->nblock; iblock++) {
    double block_energy = 0.;
    double block_energy_sqrd = 0.;
    double block_magnetization = 0.;
    double block_mag2 = 0.;
    double block_mag4 = 0.;

    for (int i = 0; i < par->L; i++)
    {
      block_spin_corrs[i] = 0.f;
    }

    for (isamp = 0; isamp < par->nsamp; isamp++) {
      int iteration = iblock*par->nsamp + isamp;
      #ifdef DRAW_GRAPHICS
      if (iteration % iteration_quotient == 0 || iteration == total_iterations-1)
      {
        draw_model(par, spin, iteration);
      }
      #endif

      accept += update(par, spin, iteration == 0 ? 1 : 0);
      
      measure(par, spin, &sample_energy, &sample_magnetization);
    }

    write_config(par, spin, fname);
    result_t r = result(par, block_energy, block_energy_sqrd, block_magnetization, block_mag2, block_mag4, par->nsamp, 0);
    datafile_write_block_results(data_file, r, iblock);
  }
  result(par, energy, energy_sqrd, magnetization, mag2, mag4, par->nblock * par->nsamp, 1);

  fclose(data_file);

  acc = accept * 100.0 / (L2 * par->nblock * par->nsamp);
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

  init_tables(par);

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
    walk_buff = realloc(walk_buff, par->N * sizeof(int));
    par->ntherm = 1000;

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
    printf("Usage: %s L=16 T=2.26\n", argv[0]);
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
