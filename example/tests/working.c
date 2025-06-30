#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct { // mpi data from simple3.c
  sc_MPI_Comm         mpicomm;
  int                 mpirank;
  int                 mpisize;
} mpi_context_t;

void debug_info(FILE* file, p8est_t *p8est, p8est_mesh_t *mesh, mpi_context_t *mpi){
  FILE * file2_debug = fopen("debug.txt", "a");

    //debug information for octants and ghost neighbors
    fprintf(file, "Rank %d: local_num_quadrants = %d\n", mpi->mpirank, mesh->local_num_quadrants);
    fprintf(file, "Rank %d: ghost_num_quadrants = %d\n", mpi->mpirank, mesh->ghost_num_quadrants);
    fprintf(file, "Rank %d: global_num_quadrants = %lld\n", mpi->mpirank, p8est->global_num_quadrants);
    fprintf(file, "Rank %d: mpisize = %d\n", mpi->mpirank, mpi->mpisize);
    fprintf(file, "\n");

    // Print all quad_to_quad values once
    fprintf(file2_debug, "Rank %d: quad_to_quad array size = %d\n", mpi->mpirank, P8EST_FACES * mesh->local_num_quadrants);
    for(int j = 0; j < P8EST_FACES * mesh->local_num_quadrants; ++j) {
        fprintf(file2_debug, "Rank %d: mesh->quad_to_quad[%d] = %d\n", mpi->mpirank, j, mesh->quad_to_quad[j]);
    }
    fprintf(file2_debug, "\n");
    fclose(file2_debug);
}

int main(int argc, char **argv) {
    int mpiret;
    mpi_context_t mpi_context, *mpi = &mpi_context;
    p8est_t *p8est;
    p8est_connectivity_t *connectivity;
    int refine_level, nX, nY, nZ, periodic;

    //pass refine level as command line argument

    //input args
    //cube -> c <refinment level>
    //brick -> b <refinment level> <nX> <nY> <nZ> <periodic>
    if (argc < 3) {
      return -1;
    }
    if (strcmp(argv[1],  "b") == 0){
      nX = atoi(argv[3]);
      nY = atoi(argv[4]);
      nZ = atoi(argv[5]);
      periodic = atoi(argv[6]); //1 = y 0 = n
    }
    refine_level = atoi(argv[2]);

    mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpi->mpicomm = sc_MPI_COMM_WORLD;
    mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
    SC_CHECK_MPI (mpiret);

    sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init (NULL, SC_LP_DEFAULT);

    if (strcmp(argv[1],  "b") == 0){
      connectivity = p8est_connectivity_new_brick(nX, nY, nZ, periodic, periodic, periodic);
    }
    else{
      connectivity = p8est_connectivity_new_unitcube();
    }
    
    p8est = p8est_new_ext(mpi->mpicomm, connectivity, 
                        /* Min number of quadrants */0, 
                        refine_level, 
                        /* uniform t/f */1, 
                        /* size of octant data struct */0, 
                        /* quadrant initalization function */ NULL,
                        /* pointer for global data struct */ NULL);

    P4EST_GLOBAL_INFO("Writing initial mesh to uniform_mesh_initial.vtu");
    p8est_vtk_write_file(p8est, NULL, "uniform_mesh_initial");

    //partition the mesh across the available processors
    P4EST_GLOBAL_INFO("Partitioning the mesh");
    p8est_partition(p8est, 0, NULL);

    //post partition (after the base guess from p8est_new_ext)
    P4EST_GLOBAL_INFO("Writing partitioned mesh to uniform_mesh_partitioned.vtu");
    p8est_vtk_write_file(p8est, NULL, "uniform_mesh_partitioned");

    p8est_ghost_t *ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    p8est_mesh_t *mesh = p8est_mesh_new(p8est, ghost, P8EST_CONNECT_FULL);

    //p8est_ghost_exchange_data(p8est, ghost, ghost_data); //exchange data between neighboring octants
    //ghost_data is a pointer to a data buffer


    char filename[256];
    snprintf(filename, sizeof(filename), "data%d.txt", mpi->mpirank);

    FILE *file = fopen(filename, "w");


    for (p4est_locidx_t i = 0; i < mesh->local_num_quadrants; ++i) {
        fprintf(file, "Rank %d: Octant %d (local):\n", mpi->mpirank, i);
        for (int face = 0; face < P8EST_FACES; ++face) {
            p4est_locidx_t neighbor_id = mesh->quad_to_quad[i * P8EST_FACES + face]; //size = P8EST_FACES * mesh->local_num_quadrants
            
          //condition-
          //if neighbor_id == i -> boundary since self ref for no face found
          //if neighbor_id != i && < mesh->local_num_quadrants -> locally owned octant, rank = mpi->mpirank
          //else (neighbor_id >= mesh->local_num_quadrants) -> ghost octant, rank = mesh->ghost_to_proc[diff from lnq]

            if (neighbor_id == i) {
                // Boundary face
                fprintf(file, "  Face %d: boundary face\n", face);
            } else if (neighbor_id < mesh->local_num_quadrants) {
                // Local neighbor
                fprintf(file, "  Face %d: neighbor is local octant %d (owned by rank %d)\n",
                       face, neighbor_id, mpi->mpirank);
            } else {
                // Ghost neighbor
                p4est_locidx_t ghost_index = neighbor_id - mesh->local_num_quadrants;
                int neighbor_rank = mesh->ghost_to_proc[ghost_index];
                fprintf(file, "  Face %d: neighbor is ghost octant %d (owned by rank %d)\n",
                       face, ghost_index, neighbor_rank);
            }
        }
        fprintf(file, "\n");
    }
    
    fclose(file);

    //clean up
    p8est_mesh_destroy(mesh);
    p8est_ghost_destroy(ghost);
    p8est_destroy(p8est);
    p8est_connectivity_destroy(connectivity);
    sc_finalize();
    
    sc_MPI_Barrier(mpi->mpicomm);

    mpiret = sc_MPI_Finalize();
    SC_CHECK_MPI(mpiret);

    return 0;
}