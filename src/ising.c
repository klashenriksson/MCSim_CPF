#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include "ran.h"

#include "ising.h"

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif



// Global variables
char fname[FNAMESIZE];	// Name for config files

void init_tables(Par *par)
{
  // Fix this (4). Initialize a table with exp(Delta E/T)
}

double *measure(Par *par, double *v, int *spin)
{
  // Fix this (3). Write a routine that measures energy and magnetization.
}



// void result(Par *par, ..., int ntot, int final)
//{
//  Fix this (3): Print out averages for energy, heat capacity, and
//  the absolute value of the magnetization per spin (remember to use 'fabs', not 'abs').
//  if (final)
//    printf("  --------  --------  --------\n");
//  printf(" %8f  %8f  %8f \n", energy, cv, magn);
//}


#ifndef CLU
int wrap(int x, int min, int max)
{
  if (x < min)
  {
    int diff = min-x;
    return wrap(max + diff, min, max);
  } else if (x > max)
  {
    int diff = x - max;
    return wrap(min + diff, min, max);
  }

  return x;
}

double sum_neighbor_spins(Par* par, int x, int y, int my_spin, int* spin)
{
  double res = 0.0;
  const int L = par->L;
  const int xleft = x-1;
  const int xright = x+1;
  const int ydown = y-1;
  const int yup = y+1;

  for (int x_i = xleft; x_i <= xright; x_i++)
  {
    for (int y_i = ydown; y_i <= yup; y_i++)
    {
      //  skip, not a neighbor
      if (x_i == x && y_i == y)
      {
        continue;
      }

      int x_wrapped = wrap(x_i, 0, L-1);
      int y_wrapped = wrap(y_i, 0, L-1);

      int nbor_spin = spin[x_wrapped + y_wrapped * par->L];
      res += (double)my_spin*(double)nbor_spin;
    }
  }

  return res;
}

int update(Par *par, int *spin)
{
  int accept = 0;
  const int L2 = par->L * par->L;

  for (int i = 0; i < L2; i++)
  {
    const int curr_spin = spin[i];
    const int negated_spin = -curr_spin;
    const int x = i % par->L;
    const int y = i / par->L;

    double energy_before = -sum_neighbor_spins(par, x, y, curr_spin, spin);
    double energy_after = -sum_neighbor_spins(par, x, y, negated_spin, spin);
    double energy_delta = energy_after - energy_before;

    double prob = min(exp(-energy_delta/par->t), 1.0);
    double epsi = dran();

    if (epsi < prob)
    {
      spin[i] = negated_spin;
      accept = accept + 1;
    }
  }

  return accept;
}
#endif


#ifdef CLU
#ifdef TRI
#define NNN 6
#else
#define NNN 4
#endif


// Fix this (Cluster update). Write the cluster update function.

#endif





int mc(Par *par, int *spin)
{
  int i, iblock, isamp, istep, ntherm = par->ntherm;
  double t = par->t, acc, accept = 0.0, L2 = par->L * par->L;


  // *** Read in the configuration for the present parameters if already present.
  // Fix this (5): if (read_config(par, spin, fname))
  // Fix this (5):   ntherm = 0;

  // *** Write out information about the run: size, temperature,...
#ifdef CLU
  printf("2D Ising model with the Wolff cluster algorithm.\n");
#else
  printf("2D Ising model with the Metropolis algorithm.\n");
#endif

  printf("\n====    %d x %d     T = %g    ====\n", par->L, par->L, par->t);
  printf("\nntherm  nblock   nsamp   seed\n");
  printf(" %5d   %5d   %5d   %d\n", ntherm, par->nblock, par->nsamp, par->seed);

  
  printf("\n energy      cv        magn     \n");

  // Thermalize the system 
  for (i = 0; i < ntherm; i++)
    update(par, spin);

  for (iblock = 0; iblock < par->nblock; iblock++) {
    
    for (isamp = 0; isamp < par->nsamp; isamp++) {
      accept += update(par, spin);
      // Fix this (3). measure(par, ..., spin);
    }
    
    // Fix this (5): write_config(par, spin, fname);
    // Fix this (3): Results for one block. result(par, ..., par->nsamp, 0);
    
  }
  // Fix this (3). Results for the whole run. result(par, ..., par->nblock * par->nsamp, 1);

  acc = accept * 100.0 / (L2 * par->nblock * par->nsamp);
  printf("\nAcceptance: %5.2f\n", acc);
  return 1;
}



int initialize_mc(Par *par, int *spin) {
  int i, L2 = par->L * par->L;
  char *f2;

  if (!par->L) {
    printf("Give system size N!\n");
    return 0;
  }


  init_ran(par->seed);

  // Fix this (4): init_tables(par);

  sprintf(fname, "%3.3d_%5.3f", par->L, par->t);

  return mc(par, spin);
}


int read_args(Par *par, char *arg)
{
  int st;
  static int *spin = NULL;
  char *s;

  if (!strcmp(arg, "run"))
    return initialize_mc(par, spin);

  s = strchr(arg, '=');

  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'.", arg);
    return 0;
  }

  *s++ = '\0';

  if (!strcmp(arg, "L")) {
    int i, L2;
    par->L = strtod(s, NULL);



    L2 = par->L * par->L;
    spin = realloc(spin, L2 * sizeof(int));
    for (i = 0; i < L2; i++)
      spin[i] = 1;
#ifdef CLU
    par->ntherm = 1000;
#else
    par->ntherm = L2;
#endif

    return 1;
  }

  if (!strcmp(arg, "T")) {
    par->t = strtod(s, NULL);
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
    return read_config(par, spin, s);
  }


  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}



int main(int argc, char *argv[])
{
  int i, iarg;
  Par *par = malloc(sizeof(Par));

  par->L = 0;
  par->t = 2.26;
  par->nblock = 1;
#ifdef CLU
  par->nsamp = 1000;
#else
  par->nsamp = 10000;
#endif

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
