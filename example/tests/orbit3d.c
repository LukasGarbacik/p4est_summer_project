#include "orbit3d.h"

void free_particle_data(local_data_t * g){
    //frees needed on-
    //actual particle buffer in buffer struct per octant
    //buffer itself
    //bounds

    //octant id and octant data struct are allocated by p4est
}
void cleanup(local_data_t * g){
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

void populate_oct_bounds(local_data_t * g, p4est_topidx_t which_tree, p8est_quadrant_t *octant){
    octant_data_t *data = (octant_data_t *) octant->p.user_data;

    int x = octant->x;
    int y = octant->y;
    int z = octant->z;

    double vertex[3]; // x y z of lower left corner
    p8est_qcoord_to_vertex(g->connectivity, which_tree, x, y, z, vertex);

    data->bounds->x_min = vertex[0];
    data->bounds->x_max = vertex[0] + 1;

    data->bounds->y_min = vertex[1];
    data->bounds->y_max = vertex[1] + 1;

    data->bounds->z_min = vertex[2];
    data->bounds->z_max = vertex[2] + 1;
}

void populate_particles(local_data_t * g, octant_data_t * data){
    //random 3 positions within bounds, tangetized 3d velo
    sc_rand_state_t rand_object;
	for(int i = 0; i < data->buffer->count; i++){
        //particle positions
	data->buffer->particles[i].x = data->bounds->x_min + sc_rand(&rand_object);
	data->buffer->particles[i].y = data->bounds->y_min + sc_rand(&rand_object);
	data->buffer->particles[i].z = data->bounds->z_min + sc_rand(&rand_object);

	//tangential velocity in 3d
	//planet mass and position through the global pointer
	
    }
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


    data->octant_id = g->num_octants++; //saves an incrementing octant_id (per rank)

    populate_oct_bounds(g, which_tree, octant);

    populate_particles(g, data);

    print_DEBUG_particle_data(g, data);
}

void run(int argc, char **argv){

    local_data_t global, *g = &global; //rank-local data
    memset (g, 0, sizeof (*g));

    //particles per quad
    g->ppq = 1000;
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

    p8est_setup(g);
    p8est_partition(g->p8est, 0, NULL); // 0 and NULL for uniform weight distribution

    g->ghost = p8est_ghost_new(g->p8est, P8EST_CONNECT_FULL);
    g->mesh = p8est_mesh_new(g->p8est, g->ghost, P8EST_CONNECT_FULL);

    cleanup(g);
}

void print_DEBUG_particle_data(local_data_t * g, octant_data_t * data){
	for(int i = 0; i < data->buffer->count; i++){
		if(i == 0){
		       	printf("\n\n Particle data of octant:%d -- \n", data->octant_id);
			printf("Octant Bounds: x: %f-%f, y: %f-%f, z:%f-%f",
				       	data->bounds->x_min,
				       	data->bounds->x_max, 
					data->bounds->y_min,
				       	data->bounds->y_max, 
					data->bounds->z_min, 
					data->bounds->z_max);
		}
		printf("\nparticle :%d x: %.3f, y: %.3f, z: %.3f", i, data->buffer->particles[i].x, data->buffer->particles[i].y, data->buffer->particles[i].z);
		if(i == data->buffer->count - 1) printf("\n\n");
	}
}

int main(int argc, char **argv){
    run(argc, argv);
    return 0;
}
