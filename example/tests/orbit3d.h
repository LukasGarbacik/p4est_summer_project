//3d orbital header

#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_connectivity.h>
#include <p8est_geometry.h>
#include <p8est_vtk.h>
#include <stdio.h>
#include <stdint.h>
#include <sc_random.h>

typedef struct {
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
} octant_bounds_t;

typedef struct {
  sc_MPI_Comm         mpicomm;
  int                 mpirank;
  int                 mpisize;
} mpi_context_t;

typedef struct { //simple particle struct for pos/vel
    double x, y, z;
    double vx, vy, vz;
    int quadrant_id;
} particle_t;

typedef struct {
    particle_t *particles; //particle buffer per quad
    int capacity;
    int count;
} particle_buffer_t;

typedef struct {
    octant_bounds_t * bounds; //fill subdomain
    particle_buffer_t * buffer; //octant-local particles 
    int octant_id;
} octant_data_t;

//global struct holding static data (read only)
typedef struct {

    //setup data
    mpi_context_t * mpi; 
    p8est_connectivity_t * connectivity;
    p8est_t * p8est; 
    p8est_ghost_t *ghost;
    p8est_mesh_t * mesh;

    int ppq; //particles per quadrant 

    double planet_xyz[3];
    double planet_mass;
    double grav_const;
    int num_octants;
} local_data_t;


//funciton headers

void mpi_set(int argc, char **argv, local_data_t * data);
void p8est_setup(local_data_t * data);
void oct_init(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *octant);
void populate_oct_bounds(local_data_t * g, p4est_topidx_t which_tree, p8est_quadrant_t *octant);
void run(int argc, char **argv);
void free_particle_data(local_data_t * g);
void populate_particles(local_data_t * g, octant_data_t * data);
void cleanup(local_data_t * g);
void print_DEBUG_particle_data(local_data_t * g, octant_data_t * data);
