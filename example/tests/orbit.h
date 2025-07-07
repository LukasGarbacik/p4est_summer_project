typedef struct {
  sc_MPI_Comm         mpicomm;
  int                 mpirank;
  int                 mpisize;
} mpi_context_t;



mpi_context_t mpi_init(int argc, char **argv);
p4est_t * p4est_setup(mpi_context_t *mpi_context);
void cleanup(p4est_t *p4est);
void run(int argc, char **argv);



