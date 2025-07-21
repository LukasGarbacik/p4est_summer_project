
#include "orbit3d.h"


void mpi_init(int argc, char **argv, local_data_t * data) {

    mpi_context_t mpi_context;

    data->mpi = &mpi_context;

    int mpiret = sc_MPI_Init(&argc, &argv);

    data->mpimpicomm = sc_MPI_COMM_WORLD;
    sc_MPI_Comm_rank(mpi_context.mpicomm, &mpi_context.mpirank);
    sc_MPI_Comm_size(mpi_context.mpicomm, &mpi_context.mpisize);

    SC_CHECK_MPI (mpiret);

    sc_init (mpi_context.mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

    global_rand_state = (uint64_t)(time(NULL) + mpi_context.mpirank);
    global_mpi_rank = mpi_context.mpirank;

}

void p8est_setup(local_data_t * data){

    p4est_init(NULL, SC_LP_DEFAULT);

    data->connectivity = p8est_connectivity_new_brick(3, 3, 3, 1, 1, 1); //3x3x3 always periodic

    data->p8est = p8est_new_ext(mpi->mpicomm, data->connectivity, 
                        /* Min number of quadrants */0, 
                        /* refine level*/0, 
                        /* uniform t/f */1, 
                        /* size of octant data struct */sizeof(octant_data_t), 
                        /* quadrant initalization function */ oct_init,
                        /* pointer for global data struct */ data);


}

void oct_init(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *octant, void *user){
    local_data_t *g = (local_data_t *)user;

    //allocate buffer and add array in g pointer
    //populate buffer

    //random posions

    //3d orbital mechanics


}

void run(int argc, char **argv){

    local_data_t global, *g = &global; //rank-local data
    memset (g, 0, sizeof (*g));

    //particles per quad
    g->ppq = 1;

    mpi_init(argc, argv, g);

    p8est_setup(g);
    p8est_partition(g->p8est, 0, NULL); // 0 and NULL for uniform weight distribution

    g->ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    g->mesh = p8est_mesh_new(p8est, ghost, P8EST_CONNECT_FULL);


    //DEBUG messages can be put here for testing (also in oct_init)

    
}

int main(int argc, char **argv){

    run(argc, argv);
    return 0;
}