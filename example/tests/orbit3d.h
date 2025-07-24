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
    double side_length;
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

    octant_data_t ** local_octants;
    int num_octants;
} local_data_t;



//octant struct to hold dynamic data


//funciton headers

void mpi_set(int argc, char **argv, local_data_t * data);
void p8est_setup(local_data_t * data);
void oct_init(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *octant);
void pointers_init(local_data_t * data);


/*
typedef struct p4est_iter_volume_info
{
  p4est_t            *p4est;
  p4est_ghost_t      *ghost_layer;
  p4est_quadrant_t   *quad;    /**< the quadrant of the callback */
  p4est_locidx_t      quadid;  /**< id in \a quad's tree array (see
                                    p4est_tree_t) */
  p4est_topidx_t      treeid;  /**< the tree containing \a quad
}
p4est_iter_volume_info_t;

*/
