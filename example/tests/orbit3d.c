#include "orbit3d.h"

void free_particle_data(local_data_t * g){
    for(int i = 0; i < g->num_octants; i++){
        if (g->local_octants[i] && g->local_octants[i]->buffer) {
            free(g->local_octants[i]->buffer->particles);
            g->local_octants[i]->buffer->particles = NULL;
        }
        if(g->local_octants[i]){
            free(g->local_octants[i]->buffer);
            free(g->local_octants[i]->bounds);

            g->local_octants[i]->buffer = NULL;
            g->local_octants[i]->bounds = NULL;
        }
    }
    free(g->local_octants);
}

void mpi_set(int argc, char **argv, local_data_t * data) {

    int mpiret = sc_MPI_Init(&argc, &argv);

    data->mpi->mpicomm = sc_MPI_COMM_WORLD;
    sc_MPI_Comm_rank(data->mpi->mpicomm, &data->mpi->mpirank);
    sc_MPI_Comm_size(data->mpi->mpicomm, &data->mpi->mpisize);

    SC_CHECK_MPI (mpiret);

    sc_init (data->mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
}

void p8est_setup(local_data_t * data){
    p4est_init(NULL, SC_LP_DEFAULT);

    data->connectivity = p8est_connectivity_new_brick(3, 3, 3, 1, 1, 1); //3x3x3 always periodic

    data->p8est = p8est_new_ext(data->mpi->mpicomm, data->connectivity, 
                        /* Min number of quadrants */0, 
                        /* refine level*/0, 
                        /* uniform t/f */1, 
                        /* size of octant data struct */sizeof(octant_data_t), 
                        /* quadrant initalization function */ oct_init,
                        /* pointer for global data struct */ data);
    data->p8est->user_pointer = data; //allow defined callbacks to access global pointer
}

void pointers_init(local_data_t * g){
    g->local_octants = (octant_data_t **) malloc(g->p8est->local_num_quadrants * sizeof(octant_data_t *));
    for(int i = 0; i < g->mpi->mpisize; i++){
        g->local_octants[i] = NULL;
    }
}

void populate_oct_bounds(local_data_t * g, p4est_topidx_t which_tree, p8est_quadrant_t *octant){
    octant_data_t *data = (octant_data_t *) octant->p.user_data;

    int level = octant->level;
    int x = octant->x;
    int y = octant->y;
    int z = octant->z;

    double vertex[3]; // x y z of lower left corner
    p8est_qcoord_to_vertex(g->connectivity, which_tree, x, y, z, vertex);

    data->bounds->side_length = P8EST_QUADRANT_LEN(level);

    printf("\n\noct_id: %d -- ll: %d, %d, %d -- ur: %d, %d, %d\n\n", data->octant_id, vertex[0], vertex[1], vertex[2], vertex[0] + data->bounds->side_length, vertex[1] + data->bounds->side_length, vertex[2] + data->bounds->side_length);
}

void oct_init(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *octant){
    local_data_t *g = (local_data_t *) p8est->user_pointer;
    octant->p.user_data = malloc(sizeof(octant_data_t));
    octant_data_t *data = (octant_data_t *) octant->p.user_data;

    data->buffer = (particle_buffer_t *) malloc(sizeof(particle_buffer_t)); 
    data->bounds = (octant_bounds_t *) malloc(sizeof(octant_bounds_t)); 

    data->buffer->count = g->ppq;
    data->buffer->capacity = g->ppq * 2;
    data->buffer->particles = (particle_t *) malloc(g->ppq * sizeof(particle_t));


    //connects octant data to global pointer to be accessed by other octants of the same rank
    g->local_octants[g->num_octants] = data;

    data->octant_id = g->num_octants; //saves an incrementing octant_id (per rank)

    g->num_octants++;

    populate_oct_bounds(g, which_tree, octant);

    //initalize randomness
    for(int i = 0; i < data->buffer->count; i++){
        //enter the particle data here for each
        //data->buffer->particles[i]
    }
}

void run(int argc, char **argv){

    local_data_t global, *g = &global; //rank-local data
    memset (g, 0, sizeof (*g));

    //particles per quad
    g->ppq = 1;
    g->planet_xyz[0] = 0.5;
    g->planet_xyz[1] = 0.5;
    g->planet_xyz[2] = 0.5;
    g->planet_mass = 0.0167;
    g->num_octants = 0;

    mpi_context_t mpi_context; //declare on the run stack
    g->mpi = &mpi_context;
    mpi_set(argc, argv, g);

    //initalize pointers so g can access all local particles
    //this is done before oct_init, so these pointers can be accessed directly
    pointers_init(g);

    p8est_setup(g);
    p8est_partition(g->p8est, 0, NULL); // 0 and NULL for uniform weight distribution

    g->ghost = p8est_ghost_new(g->p8est, P8EST_CONNECT_FULL);
    g->mesh = p8est_mesh_new(g->p8est, g->ghost, P8EST_CONNECT_FULL);


    //DEBUG messages can be put here for testing (also in oct_init)
    free_particle_data(g);
    p8est_mesh_destroy(g->mesh);
    p8est_ghost_destroy(g->ghost);
    p8est_destroy(g->p8est);
    p8est_connectivity_destroy(g->connectivity);
    //sc_finalize();
    
    sc_MPI_Barrier(g->mpi->mpicomm);

    int mpiret = sc_MPI_Finalize();
    SC_CHECK_MPI(mpiret);
}

int main(int argc, char **argv){
    run(argc, argv);
    return 0;
}