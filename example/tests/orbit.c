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
#include <time.h>

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

    srand(time(NULL) + mpi_context.mpirank);

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

particle_t particle_single_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant){
    double lower_pos[2], upper_pos[2]; //lower left and upper right

    p4est_qcoord_t length = P4EST_QUADRANT_LEN(quadrant->level);

    // Lower-left corner
    p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x, quadrant->y, lower_pos);

    // Upper-right corner
    p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + length, quadrant->y + length, upper_pos);

    //random positioning
    double rand_x = (double)rand() / RAND_MAX;
    double rand_y = (double)rand() / RAND_MAX;

    particle_t particle;

    particle.x = lower_pos[0] + (upper_pos[0] - lower_pos[0]) * rand_x;
    particle.y = lower_pos[1] + (upper_pos[1] - lower_pos[1]) * rand_y;

    //velocity of particle (0 temp)
    particle.vx = 0;
    particle.vy = 0;

    return particle;
}

void quad_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant){
    particle_buffer_t *buf = (particle_buffer_t *) quadrant->p.user_data;
    buf->count = p_per_quad;
    buf->particles = malloc(buf->count * sizeof(particle_t));
    for(int i = 0; i < p_per_quad; i++){
        buf->particles[i] = particle_single_init(p4est, which_tree, quadrant);
    }
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