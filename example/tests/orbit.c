#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <sc_notify.h>
#include <sc_options.h>
#include <stdio.h>
#include <stdlib.h>

#include "orbit.h"

#define p_per_quad 5

mpi_context_t mpi_init(int argc, char **argv) {

    mpi_context_t mpi_context;

    mpiret = sc_MPI_Init(&argc, &argv);

    mpi_context.mpicomm = sc_MPI_COMM_WORLD;
    mpi_context.mpirank = sc_MPI_Comm_rank(mpi_context.mpicomm);
    mpi_context.mpisize = sc_MPI_Comm_size(mpi_context.mpicomm);

    SC_CHECK_MPI (mpiret);

    sc_init (mpi_context.mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

    return mpi_context;
}

p4est_t * p4est_setup(mpi_context_t *mpi_context) {
    p4est_init(NULL, SC_LP_DEFAULT);

    p4est_connectivity_t *connectivity = p4est_connectivity_new_unitsquare();

    p4est_t *p4est = p4est_new_ext(mpi_context.mpicomm, connectivity, 0, 1, 
                                    /*octant struct size*/ sizeof(particle_buffer_t), 
                                    /*quadrant init fxn */ quad_init,
                                    /*pointer for global data struct */ NULL);

    return p4est;
}

void quad_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant){
    particle_buffer_t *buf = (particle_buffer_t) quadrant->p.user_data;
    buf->count = p_per_quad;
    buf->particles = malloc(buf->count * sizeof(particle_t));
    //initalize random particles within quadrant
}


void cleanup(p4est_t *p4est) {
    p4est_destroy(p4est);
    sc_finalize();
    sc_MPI_Finalize();
}

void run(int argc, char **argv) {
    mpi_context_t mpi_context = mpi_init(argc, argv);

    p4est_t *p4est = p4est_setup(&mpi_context);

    p4est_partition(p4est, 0, NULL); // 0 and NULL for uniform weight distribution

    p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    p4est_mesh_t *mesh = p4est_mesh_new(p4est, ghost, P4EST_CONNECT_FULL)



    cleanup(p4est);

}
int main(int argc, char **argv) {
    run(argc, argv);
    return 0;
}


//main plan ->>>

// create 2d unit squre connectivity
// create mesh
// start with 4 processors for 4 squres -> then move to less processors that have more quads
// create central planet that has attractive force as particles orbit
// create random number of particels in each quadrant and fill data structures
// main loop ->>>