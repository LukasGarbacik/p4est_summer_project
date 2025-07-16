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

typedef struct { //use as prefix indexes for the prebuffer array and reallocation
    int count1[2]; //lowest rank to highest
    int count2[2]; //first slot is the count, second slot is the rank filled in to access w/ particles
    int count3[2];
} received_counts;




mpi_context_t mpi_init(int argc, char **argv);
p4est_t * p4est_setup(mpi_context_t *mpi_context);

void quad_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant);
int find_quad(particle_t * particle);

particle_t particle_single_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant);
void print_particle_positions(p4est_t * p4est, p4est_mesh_t * mesh, mpi_context_t mpi_context, FILE *file, int ln);

void run(int argc, char **argv);
void loop(p4est_t *p4est, p4est_mesh_t *mesh, mpi_context_t mpi_context, int num_steps, const char *output_dir);

void free_particles(p4est_t *p4est, p4est_mesh_t * mesh);
void cleanup(p4est_t *p4est);


//pre send
void add_particle(particle_t particle, particle_buffer_t *buffer);
void remove_particle(particle_buffer_t *buffer, int index, int ln);

void send_counts(particle_buffer_t ** data, MPI_Request *mpi_sends, mpi_context_t *mpi_context);
void receive_counts(received_counts *counts, MPI_Request *mpi_recs, mpi_context_t *mpi_context);

void reallocate_buffer(particle_buffer_t * buffer, received_counts * counts);
int get_prefix_helper(particle_buffer_t * buffer, received_counts * counts, int rank);

void receive_particles(particle_buffer_t * new_buffer, MPI_Request *mpi_recs, received_counts *counts ,  mpi_context_t * mpi_context);
void send_particles(particle_buffer_t ** send_data, MPI_Request *mpi_sends, mpi_context_t *mpi_context);

void write_particles_vtk(p4est_t *p4est, p4est_mesh_t *mesh, mpi_context_t mpi_context, int step, const char *output_dir);
//void write_particles_vtu(p4est_t *p4est, p4est_mesh_t *mesh, mpi_context_t mpi_context, int step, const char *output_dir);