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
double delta_energy_table[5];

double unscaled_delta_2_energy(int deltaE)
{
  int index = (deltaE + 8) / 4;
  if (index < 0 || index > 4)
  {
    printf("Unexpected delta energy! Can not tabulate! (%d)\n", deltaE);
    return NAN;
  }

  return delta_energy_table[index];
}

void init_tables(Par *par)
{
  double delta_Es[] = {
    -8.,
    -4.,
    0.,
    4.,
    8.,
  };

  for (int i = 0; i < sizeof(delta_Es)/sizeof(delta_Es[0]); i++)
  {
    delta_energy_table[i] = exp(-delta_Es[i]/par->t);
  }
}

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

double sum_neighbor_spins(Par* par, int x, int y, int my_spin, int* spin)
{
  double res = 0.0;
  const int L = par->L;
  const int xleft = x-1;
  const int xright = x+1;
  const int ydown = y-1;
  const int yup = y+1;

  const int neighbors[] = {
    xleft, y, 
    xright, y, 
    x, yup, 
    x, ydown
  };
  const int n_neighbors = 4;

  for (int i = 0; i < n_neighbors; i++)
  {
    const int* pos = &neighbors[i * 2];
    int x_wrapped = wrap(pos[0], 0, L-1);
    int y_wrapped = wrap(pos[1], 0, L-1);
    int nbor_spin = spin[x_wrapped + y_wrapped * par->L];
    res += (double)my_spin*(double)nbor_spin;
  }

  return res;
}

void measure(Par *par, int *spin, double* out_energy, double* out_magnetization)
{
  int L2 = par->L * par->L;

  *out_energy = 0.;
  *out_magnetization = 0.f;
  for (int i = 0; i < L2; i++)
  {
    int x = i % par->L;
    int y = i / par->L;
    double energy_per_spin = 0.5 * -sum_neighbor_spins(par, x, y, spin[i], spin); //Since we double count we divide by half
    *out_energy += energy_per_spin;
    *out_magnetization += spin[i];
  }
}



void result(Par *par, double energy, double energy_sqrd, double magnetization, int ntot, int final)
{
  int L2 = par->L * par->L;
  double ntot_d = (double)ntot;
  double l2_ntot_d = (double)(L2 * ntot);
  double e = energy/l2_ntot_d;
  double Cv = (1.0 / (par->t * par->t)) * (energy_sqrd/ntot_d - energy*energy/(ntot_d*ntot_d))/(double)L2;
  double m = magnetization/l2_ntot_d;

  if (final)
  {
    printf("  --------  --------  --------\n");

    //also write results to file.
    char resname[256] = {0};
    sprintf(resname, "results_T_%8f_L_%d.txt", par->t, par->L);
    FILE* res_file = fopen(resname, "a+");
    if (res_file)
    {
      fprintf(res_file, "%8f %d %8f %8f %8f\n", par->t, par->L, e, Cv, m);
      fclose(res_file);
    }
  }
  printf(" %8f  %8f  %8f \n", e, Cv, m);
}


#ifndef CLU
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

    double prob = min(unscaled_delta_2_energy(energy_delta), 1.0);
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
  double energy = 0.0, energy_sqrd = 0.0, magnetization = 0.0;


  // *** Read in the configuration for the present parameters if already present.
  if (read_config(par, spin, fname))
    ntherm = 0;

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
    double block_energy = 0.;
    double block_energy_sqrd = 0.;
    double block_magnetization = 0.;

    for (isamp = 0; isamp < par->nsamp; isamp++) {
      accept += update(par, spin);
      // Fix this (3). measure(par, ..., spin);
      
      double sample_energy = 0.0, sample_magnetization = 0.0;
      measure(par, spin, &sample_energy, &sample_magnetization);

      block_energy += sample_energy;
      block_energy_sqrd += sample_energy * sample_energy;
      block_magnetization += fabs(sample_magnetization); //We only care about the abolute magnetization, not direction
    }

    energy += block_energy;
    energy_sqrd += block_energy_sqrd;
    magnetization += block_magnetization;

    write_config(par, spin, fname);
    result(par, block_energy, block_energy_sqrd, block_magnetization, par->nsamp, 0);
    
  }
  // Fix this (3). Results for the whole run. result(par, ..., par->nblock * par->nsamp, 1);
  result(par, energy, energy_sqrd, magnetization, par->nblock * par->nsamp, 1);

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

  init_tables(par);

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
