//3d orbital header

#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_connectivity.h>
#include <p8est_mesh.h>
#include <p8est_ghost.h>
#include <p8est_geometry.h>
#include <p8est_vtk.h>
#include <sc.h>
#include <sc_random.h>
#include <sys/stat.h>
#include <mpi.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

#define total_octants 27
#define global_bound 3.0

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
    int octant_id;
} particle_t;

typedef struct {
    particle_t *particles; //particle buffer per quad
    int capacity;
    int count;
} particle_buffer_t;

typedef struct {
    octant_bounds_t bounds; //fill subdomain
    particle_buffer_t * buffer; //octant-local particles 
    int octant_id;
    int mpirank; //makes rank public to other octants when inserting ghost particles (corner funciton)
} octant_data_t;

//global struct holding static data (read only)
typedef struct {

    //mpi and p4est objects
    mpi_context_t * mpi; 
    p8est_connectivity_t * connectivity;
    p8est_t * p8est; 
    p8est_ghost_t *ghost;
    p8est_mesh_t * mesh;
    sc_array_t * ghost_data;

    //static information
    int ppq; //particles per quadrant 
    double planet_xyz[3];
    double planet_mass;
    double grav_const;
    int num_octants;
    int num_steps;
    double timestep;
    int total_considered_oct;
    int no; //ghost estimate for number of octants (close upper bound assuming uniform octant distribution)

    //send/recv
    octant_data_t * outgoing_local; //local array of octant data to be loaded
    octant_data_t ** outgoing_ghost; // outgoing_ghost[i] (i) is rank of recv, outgoing_ghost[i]->octants[j].octant_id is octant_id of recv octant over ghost

    int **send_counts;
    int **recv_counts;
} local_data_t;

//setup
void mpi_set(int argc, char **argv, local_data_t * data);
void p8est_setup(local_data_t * data);
void oct_init(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *octant);
//main functionality
void populate_oct_bounds(local_data_t * g, p4est_topidx_t which_tree, p8est_quadrant_t *octant);
void populate_particles(local_data_t * g, octant_data_t * data);
void populate_send_buffers(local_data_t * g);
void insert_particle_into_outgoing(local_data_t * g, p8est_quadrant_t *oct, particle_t * p, int new_id);
void insert_particle(particle_t * p, particle_buffer_t * buffer);
void insert_ghost_particle(local_data_t * g, p8est_quadrant_t *oct, particle_t * p);
void remove_particle(particle_buffer_t * buffer, int index);
void combine_local(local_data_t * g);
void do_dynamics(local_data_t * g);
void send_ghost(local_data_t * g);
void send_recv_counts(local_data_t * g);
void send_recv_particles(local_data_t * g);
//helper
bool in_bounds(particle_t * p, octant_bounds_t * bounds);
int query_oct_id(local_data_t * g, particle_t * p);
//main run
void run(int argc, char **argv);
void loop(local_data_t * g, const char *output_dir);
//cleanup
void cleanup(local_data_t * g);
void free_particle_data(local_data_t * g);
void free_local(local_data_t * g);
void free_ghost(local_data_t * g);
//output
//void print_DEBUG_particle_data(const octant_data_t * data);
//void print_DEBUG_send_data(const local_data_t * g);
void write_vtk(local_data_t * g, const char *output_dir, int cur_step);
void print_long_DEBUG(local_data_t * g, int loop_num);
void print_prefix_matrix(local_data_t *g);
