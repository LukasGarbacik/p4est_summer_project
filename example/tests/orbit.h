typedef struct {
  sc_MPI_Comm         mpicomm;
  int                 mpirank;
  int                 mpisize;
} mpi_context_t;

typedef struct { //simple particle struct for pos/vel
    double x, y;
    double vx, vy;
    int quadrant_id; //same value as mpi rank
} particle_t;

typedef struct {
    particle_t *particles; //particle buffer per quad
    int capacity; //change to realloc (dynamic resizing) later
    int count; // number of active owned particles
} particle_buffer_t;





mpi_context_t mpi_init(int argc, char **argv);
p4est_t * p4est_setup(mpi_context_t *mpi_context);

void quad_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant);
int find_quad(particle_t * particle);

particle_t particle_single_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant);
void print_particle_positions(p4est_t * p4est, p4est_mesh_t * mesh, mpi_context_t mpi_context, FILE *file);

void run(int argc, char **argv);
void loop(p4est_t *p4est, p4est_mesh_t *mesh, mpi_context_t mpi_context, int num_steps, FILE *file);

void free_particles(p4est_t *p4est, p4est_mesh_t * mesh);
void cleanup(p4est_t *p4est);


//pre send
void add_particle(particle_t particle, particle_buffer_t *buffer);
void remove_particle(particle_buffer_t *buffer, int index);


