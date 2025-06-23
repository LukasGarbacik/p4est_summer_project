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

int main(int argc, char **argv) {
    int mpiret;
    mpi_context_t mpi_context, *mpi = &mpi_context;
    p8est_t *p8est;
    p8est_connectivity_t *connectivity;
    int                 refine_level;

    /* Check for command line argument */
    if (argc < 2) {
      printf("Usage: %s <refine_level>\n", argv[0]);
      return -1;
    }
    refine_level = atoi(argv[1]);

    mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpi->mpicomm = sc_MPI_COMM_WORLD;
    mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
    SC_CHECK_MPI (mpiret);

    sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init (NULL, SC_LP_DEFAULT);

    /* Create a connectivity that is a single unit cube */
    connectivity = p8est_connectivity_new_unitcube();
    
    /* Create the p8est forest.
     * The `uniform` flag (the 5th argument, set to 1) tells p4est to
     * automatically create a uniformly refined mesh up to `min_level`
     * (the 4th argument). This is exactly what you want. */
    p8est = p8est_new_ext(mpi->mpicomm, connectivity,
                          0,            /* min_quadrants */
                          refine_level, /* min_level */
                          1,            /* uniform */
                          0,            /* user_data_size */
                          NULL,         /* init_fn */
                          NULL);        /* user_pointer */

    /* Save the initial mesh to a file */
    P4EST_GLOBAL_INFO("Writing initial mesh to uniform_mesh_initial.vtu");
    p8est_vtk_write_file(p8est, NULL, "uniform_mesh_initial");

    /* Partition the mesh across the available processors */
    P4EST_GLOBAL_INFO("Partitioning the mesh");
    p8est_partition(p8est, 0, NULL);

    /* Save the partitioned mesh to a file */
    P4EST_GLOBAL_INFO("Writing partitioned mesh to uniform_mesh_partitioned.vtu");
    p8est_vtk_write_file(p8est, NULL, "uniform_mesh_partitioned");

    /* Clean up resources */
    p8est_destroy(p8est);
    p8est_connectivity_destroy(connectivity);
    sc_finalize();

    mpiret = sc_MPI_Finalize();
    SC_CHECK_MPI(mpiret);

    return 0;
}