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
                                    /*octant struct size*/0, 
                                    /*quadrant init fxn */0,
                                    /*pointer for global data struct */ NULL);

    return p4est;
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