#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#include "ran.h"
#include "percolation.h"
#include "int_queue.h"

#include "data.h"

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

// Global variables
char fname[FNAMESIZE];	// Name for config files

typedef struct lattice_site {
  int flags;
  int vis_x,vis_y,vis_z;
} lattice_site_t;

const int nearest_neighbors[] = {
  #if D == 2
  1,0,
  -1,0,
  0,1,
  0,-1,
  #elif D == 3
  1,0,0,
  -1,0,0,
  0,1,0,
  0,-1,0,
  0,0,1,
  0,0,-1,
  #endif
};
#if D == 2
const int n_neighbors = 4;

#elif D == 3
const int n_neighbors = 6;
#endif

const int F_OCCUPIED = 1 << 0;
const int F_VISITED = 1 << 1;
const int Y_AXIS = 0;
const int X_AXIS = 1;
const int Z_AXIS = 2;

void get_neighbor_xyz(int n_idx, int* x, int* y, int* z)
{
  #if D == 2
  *x = nearest_neighbors[n_idx*2];
  *y = nearest_neighbors[n_idx*2 + 1];
  *z = 0;
  #elif D == 3
  *x = nearest_neighbors[n_idx * 3];
  *y = nearest_neighbors[n_idx * 3 + 1];
  *z = nearest_neighbors[n_idx * 3 + 2];
  #endif
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

int has_flags(Par* par, lattice_site_t* lattice, int x, int y, int z, int f)
{
  x = wrap(x, 0, par->L-1);
  y = wrap(y, 0, par->L-1);
  z = wrap(z, 0, par->L-1);
  return lattice[x + y * par->L + z * par->L * par->L].flags & f;
}

void enable_flags(Par* par, lattice_site_t* lattice, int x, int y, int z, int f)
{
  x = wrap(x, 0, par->L-1);
  y = wrap(y, 0, par->L-1);
  z = wrap(z, 0, par->L-1);
  lattice[x + y * par->L + z * par->L * par->L].flags |= f;
}

void disable_flags(Par* par, lattice_site_t* lattice, int x, int y, int z, int f)
{
  x = wrap(x, 0, par->L-1);
  y = wrap(y, 0, par->L-1);
  z = wrap(z, 0, par->L-1);
  lattice[x + y * par->L + z * par->L * par->L].flags &= ~f;
}

void get_lattice_vis_pos(Par* par, lattice_site_t* lattice, int x,int y,int z, int* out_x, int* out_y, int* out_z)
{
  x = wrap(x, 0, par->L-1);
  y = wrap(y, 0, par->L-1);
  z = wrap(z, 0, par->L-1);
  lattice_site_t site = lattice[x + y * par->L + z * par->L * par->L];
  *out_x = site.vis_x;
  *out_y = site.vis_y;
  *out_z = site.vis_z;
}

void set_lattice_vis_pos(Par* par, lattice_site_t* lattice, int x,int y,int z, int vis_x, int vis_y, int vis_z)
{
  x = wrap(x, 0, par->L-1);
  y = wrap(y, 0, par->L-1);
  z = wrap(z, 0, par->L-1);
  lattice_site_t* site = &lattice[x + y * par->L + z * par->L * par->L];
  site->vis_x = vis_x;
  site->vis_y = vis_y;
  site->vis_z = vis_z;
}

void idx2xyz(Par* par, int i, int* x, int* y, int* z)
{
  #if D == 3
  *z = i / (par->L * par->L);
  i -= *z * par->L * par->L;
  *y = i / par->L;
  *x = i % par->L;
  #elif D == 2
  *x = i % par->L;
  *y = i / par->L;
  *z = 0;
  #endif
}

int xyz2idx(Par* par, int x, int y, int z)
{
  x = wrap(x, 0, par->L-1);
  y = wrap(y, 0, par->L-1);
  z = wrap(z, 0, par->L-1);
  #if D == 3
  return x + y*par->L + z * par->L * par->L;
  #elif D == 2
  return x + y*par->L;
  #endif
}

int inside(Par* par, int x)
{
  return (x >= 0 && x < par->L);
}

void print_lattie(Par* par, lattice_site_t* lattice)
{
  printf("--------------\n");
  #if D == 2
  for (int y = 0; y < par->L; y++)
  {
    for(int x = 0; x < par->L; x++)
    {
      if(has_flags(par, lattice, x,y,0, F_OCCUPIED))
      {
        printf("x");
      } else {
        printf("*");
      }
    }
    printf("\n");
  }
  #else
  printf("not supported for requested dimension\n");
  #endif
}

int check_percolation(Par* par, lattice_site_t* lattice, int_queue_t* queue, int start_x, int start_y, int start_z, int dir)
{
  if (has_flags(par, lattice, start_x, start_y, start_z, F_VISITED) || !has_flags(par, lattice, start_x, start_y, start_z, F_OCCUPIED))
  {
    return 0;
  }

  enable_flags(par, lattice, start_x,start_y,start_z, F_VISITED);
  set_lattice_vis_pos(par, lattice, start_x, start_y, start_z, start_x, start_y, start_z);
  int idx = xyz2idx(par, start_x, start_y, start_z);

  int_queue_push_back(queue, idx);

  while (int_queue_size(queue) > 0)
  {
    idx = int_queue_pop_front(queue);
    int x,y,z; //guaranteed to be within [0,L-1]
    idx2xyz(par, idx, &x, &y, &z);

    int curr_stored_x, curr_stored_y, curr_stored_z;
    get_lattice_vis_pos(par, lattice, x,y,z, &curr_stored_x, &curr_stored_y, &curr_stored_z);

    //check each nearest neighbor
    for (int i = 0; i < n_neighbors; i++)
    {
      int nx, ny, nz; //relative coords
      get_neighbor_xyz(i, &nx, &ny, &nz);

      int nbor_x = x + nx; //guaranteed to be within [-1,L]
      int nbor_y = y + ny;
      int nbor_z = z + nz;

      int nx_rel = curr_stored_x + nx;
      int ny_rel = curr_stored_y + ny;
      int nz_rel = curr_stored_z + nz;

      // If not occupied, it is not of interest
      if (!has_flags(par, lattice, nbor_x, nbor_y, nbor_z, F_OCCUPIED))
      {
        continue;
      }

      if (has_flags(par, lattice, nbor_x, nbor_y, nbor_z, F_VISITED))
      {
        // We might have a percolation here!
        int stored_x, stored_y, stored_z;
        get_lattice_vis_pos(par, lattice, nbor_x, nbor_y, nbor_z, &stored_x, &stored_y, &stored_z);
        int dx = stored_x - nx_rel;
        int dy = stored_y - ny_rel;
        int dz = stored_z - nz_rel;

        if (abs(dx) == par->L || abs(dy) == par->L || abs(dz) == par->L)
        {
          print_lattie(par, lattice);
          return 1;
        }
      }
      else
      {
        // No percolation
        set_lattice_vis_pos(par, lattice, nbor_x, nbor_y, nbor_z, nx_rel, ny_rel, nz_rel);
        enable_flags(par, lattice, nbor_x, nbor_y, nbor_z, F_VISITED);
        int_queue_push_back(queue, xyz2idx(par, nbor_x, nbor_y, nbor_z));
      }
    }

  }

  return 0;
}

void measure(Par *par, lattice_site_t* lattice, int_queue_t* queue, int* out_did_percolate)
{
  *out_did_percolate = 0;

  int y_q = 0;
  int z_loop = 0;

  for (int x_loop = 0; x_loop < par->L; x_loop++)
  {
    #if D == 3
    for (z_loop = 0; z_loop < par->L; z_loop++)
    {
    #endif

    if (check_percolation(par, lattice, queue, x_loop, y_q, z_loop, Y_AXIS))
    {
      *out_did_percolate = 1;
      return;
    }

    /*if (has_flags(par, lattice, x_loop, y_q, z_loop, F_OCCUPIED) && !has_flags(par, lattice, x_loop, y_q, z_loop, F_VISITED))
    {
      int idx = xyz2idx(par, x_loop, y_q, z_loop);

      enable_flags(par, lattice, x_loop,y_q,z_loop, F_VISITED);
      set_lattice_vis_pos(par, lattice, x_loop, y_q, z_loop, x_loop, y_q, z_loop);

      int_queue_push_back(queue, idx);

      while(int_queue_size(queue) > 0)
      {
        int idx = int_queue_pop_front(queue);
        int x,y,z; //actual coordinate on lattice
        idx2xyz(par, idx, &x,&y,&z);

        int vis_x, vis_y, vis_z; //coordinate stored at the site
        get_lattice_vis_pos(par, lattice, x,y,z, &vis_x, &vis_y, &vis_z);

        for (int i = 0; i < n_neighbors; i++)
        {
          int nx,ny,nz;
          get_neighbor_xyz(i, &nx, &ny, &nz);

          //  neighbour coords relative to position stored at site
          int vrel_nx = vis_x + nx;
          int vrel_ny = vis_y + ny;
          int vrel_nz = vis_z + nz;

          //  neighbour coords relative to actual coordinate
          nx += x;
          ny += y;
          nz += z;

          if (has_flags(par, lattice, nx,ny,nz, F_OCCUPIED))
          {
            if (has_flags(par, lattice, nx,ny,nz, F_VISITED))
            {
              //  get stored coordinates on the visited lattice
              int nv_x, nv_y, nv_z;
              get_lattice_vis_pos(par, lattice, nx, ny, nz, &nv_x, &nv_y, &nv_z);
              //if diff in y between neighbors position rel. to current stored coord and neighbors stored y coord is L
              if (abs(vrel_ny - nv_y) == par->L)
              {
                *out_did_percolate = 1;
                return;
              }
            }
            else
            {
              enable_flags(par, lattice, nx,ny,nz, F_VISITED);
              int_queue_push_back(queue, xyz2idx(par,nx,ny,nz));

              set_lattice_vis_pos(par, lattice, nx,ny,nz, vrel_nx, vrel_ny, vrel_nz);
            }
          }
        }
      }
    }*/

    #if D == 3
    }
    #endif
  }
}

result_t result(Par *par, double percs, int ntot)
{
  double perc_prob = percs / ntot;
  printf("  %8f\n", perc_prob);

  result_t r;
  r.perc_prob = perc_prob;

  return r;
}

int update(Par* par, lattice_site_t* lattice)
{
  #if D == 2
  int Ls = par->L * par->L;
  #elif D == 3
  int Ls = par->L * par->L * par->L;
  #endif
  for (int i = 0; i < Ls; i++)
  {
    double r = dran();
    int x,y,z;
    idx2xyz(par, i, &x,&y,&z);
    disable_flags(par, lattice, x,y,z, F_VISITED | F_OCCUPIED);
    
    if (r <= par->p)
      enable_flags(par, lattice, x,y,z, F_OCCUPIED);
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

int mc(Par *par, lattice_site_t *lattice)
{
  int i, iblock, isamp, istep, ntherm = par->ntherm;

  char datafilename[256] = {0};
  datafile_get_filename(par, datafilename);
  FILE* data_file = fopen(datafilename, "a");
  if(!data_file)
  {
    fprintf(stderr, "Error! Expected datafile to be created. Aborting.\n");
    return -1;
  }

  // *** Write out information about the run: size, temperature,...
  #if D == 2
  printf("2D Lattice Percolation model.\n");
  #elif D == 3
  printf("3D Lattice Percolation model.\n");
  #endif

  printf("\n====    L = %d  p = %8f  ====\n", par->L, par->p);
  printf("\nntherm  nblock   nsamp   seed\n");
  printf(" %5d   %5d   %5d   %d\n", ntherm, par->nblock, par->nsamp, par->seed);
  printf("\n  Perc. prob\n");

  #if D == 2  
  int Ls = par->L * par->L;
  #elif D == 3
  int Ls = par->L * par->L * par->L;
  #endif
  int_queue_t measure_queue = int_queue_create(Ls*4); //idx, vis coord

  for (iblock = 0; iblock < par->nblock; iblock++) {
    double block_percs = 0;
    for (isamp = 0; isamp < par->nsamp; isamp++) {
      double sample_s2 = 0.f;
      update(par, lattice);
      
      int did_perc = 0;
      int_queue_empty(&measure_queue);
      measure(par, lattice, &measure_queue, &did_perc);

      block_percs += did_perc;
    }

    result_t r = result(par, block_percs, par->nsamp);
    datafile_write_block_results(data_file, r, iblock);
  }

  fclose(data_file);
  int_queue_destroy(&measure_queue);

  return 1;
}

int initialize_mc(Par *par, lattice_site_t *lattice) {
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
  static lattice_site_t *lattice = NULL;
  char *s;

  if (!strcmp(arg, "run"))
    return initialize_mc(par, lattice);

  s = strchr(arg, '=');

  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'.", arg);
    return 0;
  }

  *s++ = '\0';

  if (!strcmp(arg, "L")) {
    par->L = strtod(s, NULL);
    #if D == 2
    int Ls = par->L * par->L;
    #elif D == 3
    int Ls = par->L * par->L * par->L;
    #endif
    lattice = realloc(lattice, Ls * sizeof(lattice_site_t));
    par->ntherm = 0;

    return 1;
  }

  if (!strcmp(arg, "p")) {
    par->p = strtod(s, NULL);
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

  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}

int main(int argc, char *argv[])
{
  int i, iarg;
  Par *par = malloc(sizeof(Par));

  par->L = 0;
  par->p = 0.5;
  par->nblock = 1;
  par->nsamp = 10000;

  if (argc == 1) {
    printf("Usage: %s L=16\n", argv[0]);
    printf("Optional arguments (with defaults) nblock=%d nsamp=%d ntherm=%d seed=%d, p=%lf\n",
	   par->nblock, par->nsamp, par->ntherm, par->seed, 0.5);
    exit(EXIT_SUCCESS);
  }

  // Interpret the commands given in the argument list.
  for (iarg = 1; iarg < argc; iarg++)
    if (!read_args(par, argv[iarg]))
      exit(EXIT_FAILURE);

  free(par);
  return 0;
}
