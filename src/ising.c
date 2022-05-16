#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#include "ran.h"
#include "ising.h"
#include "data.h"

#ifdef CORR
#include "corr.h"
#endif

#ifdef DRAW_GRAPHICS
#include <g2.h>
#include <g2_gd.h>
#endif

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
    x, ydown,
    #ifdef TRI
    xleft, yup,
    xright, ydown,
    #endif
  };

  for (int i = 0; i < sizeof(neighbors)/sizeof(int)/2; i++)
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


result_t result(Par *par, double energy, double energy_sqrd, double magnetization, double mag2, double mag4, int ntot, int final)
{
  double L2 = (double)par->L * (double)par->L;
  double ntot_d = (double)ntot;

  double E = energy/ntot_d;
  double E2 = energy_sqrd/ntot_d;

  double e = E/L2;

  double Cv = (1.0 / (par->t * par->t)) * (E2 - E*E) / L2;

  double M = magnetization/ntot_d;
  double M2 = mag2/ntot_d;
  double M4 = mag4/ntot_d;

  double m = M/L2;
  double m2 = M2/(L2*L2);
  double m4 = M4/(L2*L2*L2*L2);

  if (final)
  {
    printf("  --------  --------  --------\n");
  }
  printf(" %8f  %8f  %8f %8f\n", e, Cv, m, m2);

  result_t r;
  r.e = e;
  r.c = Cv;
  r.m = m;
  r.m2 = m2;
  r.m4 = m4;

  return r;
}

#ifdef DRAW_GRAPHICS
void draw_model(Par* par, int* spins, int iteration)
{
  char filename[256] = {0};
  sprintf(filename, "model_%f_%d_%d.png", par->t, par->L, iteration);


  const int image_width = 800;
  const int image_height = 800;
  const float spin_width = (float)image_width/(float)par->L;
  const float spin_height = (float)image_height/(float)par->L;
  const float spin_r = sqrt(0.25 * spin_width * spin_width + 0.25 * spin_height + spin_height);

  int surf = g2_open_gd(filename, image_width, image_height, g2_gd_png);

  const int L2 = par->L * par->L;
  for (int i = 0; i < L2; i++)
  {
    const float x = (float)(i % par->L) * spin_width;
    const float y = (float)(i / par->L) * spin_height;

    int spin = spins[i];
    int color = spin + 2;

    g2_pen(surf, color);
    g2_filled_circle(surf, x + spin_r, y + spin_r, spin_r);
  }

  g2_close(surf);
}
#endif

#ifndef CLU
int update(Par *par, int *spin, int draw)
{
  int accept = 0;
  const int L2 = par->L * par->L;

  #ifdef DRAW_GRAPHICS
  if (draw)
    draw_model(par, spin, 0);
  #endif

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
      #ifdef DRAW_GRAPHICS
      if (draw)
        draw_model(par, spin, i+1);
      #endif
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


int update(Par* par, int* spin, int draw)
{
  #ifdef DRAW_GRAPHICS
  if (draw)
  {
    draw_model(par, spin, 0);
  }
  #endif
  int accept = 1;
  int_queue_empty(&par->queue);

  const int L2 = par->L * par->L;
  const int i_start = iran() % L2;
  const int S = spin[i_start];
  spin[i_start] = -S;

  int_queue_push_back(&par->queue, i_start);

#ifdef DRAW_GRAPHICS
  int iteration = 0;
#endif
  while (int_queue_size(&par->queue) > 0)
  {
    const int i = int_queue_pop_front(&par->queue);
    const int x = i % par->L;
    const int y = i / par->L;
    const int xleft = x-1;
    const int xright = x+1;
    const int yup = y+1;
    const int ydown = y-1;

    const int neighbors[] = {
      xleft, y, 
      xright, y, 
      x, yup, 
      x, ydown,
#ifdef TRI
      xleft, yup,
      xright, ydown
#endif
    };
    const int n_neighbors = sizeof(neighbors)/sizeof(int)/2;

    for (int n_i = 0; n_i < n_neighbors; n_i++)
    {
      const int* pos = &neighbors[n_i*2];
      const int x_wrapped = wrap(pos[0], 0, par->L-1);
      const int y_wrapped = wrap(pos[1], 0, par->L-1);
      const int nbor_spin = spin[x_wrapped + y_wrapped * par->L];

      if (nbor_spin == S)
      {
        const double prob = min(1.0 - exp(-2.0/par->t), 1.0);
        const double epsi = dran();
        if (epsi < prob)
        {
          int nbor_idx = x_wrapped + y_wrapped * par->L;
          spin[nbor_idx] = -spin[nbor_idx];
          int_queue_push_back(&par->queue, nbor_idx);
          accept += 1;
          #ifdef DRAW_GRAPHICS
          if (draw)
            draw_model(par, spin, iteration++);
          #endif     
        }
      }
    }
  }

  return accept;
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
      
      double sample_energy = 0.0, sample_magnetization = 0.0;
      measure(par, spin, &sample_energy, &sample_magnetization);

      #ifdef CORR
      corr_step(&corr, sample_energy);
      #endif

      // spin correlation
      for (int x = 0; x < par->L; x++)
      {
        double spin_corr = 0.;
        for (int x_sweep = 0; x_sweep < par->L; x_sweep++)
        {
          int y = uran() % par->L;
          int idx1 = x_sweep + y*par->L;
          int idx2 = wrap(x_sweep+x, 0, par->L-1) + y*par->L;
          int spin_1 = spin[idx1];
          int spin_2 = spin[idx2];
          spin_corr += spin_1*spin_2;
        }

        spin_corr /= (par->L * par->nsamp);
        block_spin_corrs[x] += spin_corr;
      }

      block_energy += sample_energy;
      block_energy_sqrd += sample_energy * sample_energy;
      block_magnetization += fabs(sample_magnetization); //We only care about the abolute magnetization, not direction
      block_mag2 += sample_magnetization * sample_magnetization;
      block_mag4 += sample_magnetization * sample_magnetization * sample_magnetization * sample_magnetization;
    }

    energy += block_energy;
    energy_sqrd += block_energy_sqrd;
    magnetization += block_magnetization;
    mag2 += block_mag2;
    mag4 += block_mag4;

    write_config(par, spin, fname);
    result_t r = result(par, block_energy, block_energy_sqrd, block_magnetization, block_mag2, block_mag4, par->nsamp, 0);
    datafile_write_block_results(data_file, r, iblock);
    datafile_write_block_spin_correlations(data_file, block_spin_corrs, par->L);
  }
  result(par, energy, energy_sqrd, magnetization, mag2, mag4, par->nblock * par->nsamp, 1);

  fclose(data_file);

  acc = accept * 100.0 / (L2 * par->nblock * par->nsamp);
  printf("\nAcceptance: %5.2f\n", acc);

  #ifdef CORR
  create_if_not_exists("correlations/");
  char corrname[256] = {0};
  sprintf(corrname, "correlations/correlations_T_%8f_L_%d.txt", par->t, par->L);
  FILE* corr_file = fopen(corrname, "w");
  if (corr_file)
  {
    double correlations[200] = {0};
    corr_compute(&corr, correlations, 200);

    for(int i = 0; i < 200; i++)
    {
      fprintf(corr_file, "%8f %d %d %5.2f\n", par->t, par->L, i, correlations[i]);
    }

    fclose(corr_file);
  }
  #endif

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

  check_data_file(par);

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
