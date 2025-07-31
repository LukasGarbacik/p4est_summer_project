#include "orbit3d.h"

static sc_rand_state_t rand_object;

void free_particle_data(local_data_t * g){
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        if(data){
            if(data->buffer){
                if(data->buffer->particles) free(data->buffer->particles);
                free(data->buffer);
            }
            if(data->bounds) free(data->bounds);
            oct->p.user_data = NULL;
        }
    }
}

void cleanup(local_data_t * g){
    //free_particle_data(g); having an issue here
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
    rand_object = (uint64_t)(time(NULL) + g->mpi->mpirank);
    double dx, dy, dz;
	for(int i = 0; i < data->buffer->count; i++){
        //particle positions
        particle_t * p = &data->buffer->particles[i];
	    p->x = data->bounds->x_min + sc_rand(&rand_object);
        printf("RAND DEBUG EXAMPLE: %.4f", sc_rand(&rand_object));
        printf("\nx bound:%.2f-%.2f p->x: %.3f\n", data->bounds->x_min, data->bounds->x_max, p->x);
	    p->y = data->bounds->y_min + sc_rand(&rand_object);
        printf("y bound:%.2f-%.2f p->y: %.3f\n", data->bounds->y_min, data->bounds->y_max, p->y);
	    p->z = data->bounds->z_min + sc_rand(&rand_object);
        printf("z bound:%.2f-%.2f p->z: %.3f\n", data->bounds->z_min, data->bounds->z_max, p->z);

	    //tangential velocity in 3d
        dx = p->x - g->planet_xyz[0];
        dy = p->y - g->planet_xyz[1];
        dz = p->z - g->planet_xyz[2];
        double r = sqrt(dx*dx + dy*dy + dz*dz); //dist from planet
        double v = sqrt(g->grav_const * g->planet_mass / r); //magnitude of velocity
        double ax = 0.7, ay = 1, az = 0.4; //arbitrary vector 
        if (fabs(ax - dx) < 1e-6 && fabs(ay - dy) < 1e-6) {
            //if arb-vec is parallel w/ radial -> pick new arb-vec
            ax = 0.2; ay = 0.5; az = 1;
        }
           // tangent in direction (cross product of radial vector and arbitrary vector)
           // 3 x 3 matrix cross product equations
            double x_vec = dy * az - dz * ay;
            double y_vec = dz * ax - dx * az;
            double z_vec = dx * ay - dy * ax;
            double mag = sqrt(x_vec*x_vec + y_vec*y_vec + z_vec*z_vec);

            p->vx = v * x_vec / mag;
            p->vy = v * y_vec / mag;
            p->vz = v * z_vec / mag;

            if(!in_bounds(p, data->bounds)) printf("INIT ISSUE WITH BOUNDS");
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
}


bool in_bounds(particle_t * p, octant_bounds_t * bounds){
    if(
        p->x >= bounds->x_min && p->x <= bounds->x_max &&
        p->y >= bounds->y_min && p->y <= bounds->y_max &&
        p->z >= bounds->z_min && p->z <= bounds->z_max
    ) return true;

    return false;
}
int query_oct_id(local_data_t * g, particle_t * p){
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        octant_data_t * data = (octant_data_t *) oct->p.user_data;

        if(in_bounds(p, data->bounds)){
            printf("\nnew oct: %d, particle (%.2f, %.2f, %.2f) being sent", data->octant_id, p->x, p->y, p->z);
            printf("\nbounds- x: %.2f, y:%.2f, z:%.2f \n", data->bounds->x_min,data->bounds->y_min, data->bounds->z_min);
            return data->octant_id;
        }
    }
    return -1; //send to ghost
}

void insert_particle(particle_t * p, particle_buffer_t * buffer){
    if(buffer == NULL || buffer->particles == NULL){
        return;
    }
    if(buffer->count >= buffer->capacity){
        buffer->capacity *= 2;  // Double the capacity
        buffer->particles = realloc(buffer->particles, buffer->capacity * sizeof(particle_t));
    }
    
    buffer->particles[buffer->count] = *p;
    buffer->count++;
}


void insert_particle_into_outgoing(local_data_t * g, particle_t * p, int new_id){
    if(new_id == -1){

    } //insert to ghost
    else{
        printf("inserting partice into octant: %d", new_id);
        if(!g->outgoing_local){ //allocate 2d array if not yet and init oct_id's to -2
            g->outgoing_local = (octant_data_t *) malloc(g->mesh->local_num_quadrants * sizeof(octant_data_t));
            memset(g->outgoing_local, 0, g->mesh->local_num_quadrants * sizeof(octant_data_t));
            printf("\nallocated %d bytes", g->mesh->local_num_quadrants * sizeof(octant_data_t));
            for(int i = 0; i < g->mesh->local_num_quadrants; ++i){
                g->outgoing_local[i].octant_id = -2; //dummy value
            }
        }
        printf("finished dummy set");
        //each process combines transfer particles into one buffer per octant
        for(int i = 0; i < g->mesh->local_num_quadrants; ++i){
            //wrong buffer ->
            printf("\ndoing loop");
            if(g->outgoing_local[i].octant_id != -2 && g->outgoing_local[i].octant_id != new_id) continue;
            //not the first particle in local send buffer
            if(g->outgoing_local[i].octant_id == new_id){
                printf("\nparticle about to be inserted into buffer fr");
                insert_particle(p, g->outgoing_local[i].buffer);
            }
            else{//is the first particle in local send buffer (init id)
                printf("\nDOUBLE sanity check local send");
                g->outgoing_local[i].buffer = (particle_buffer_t *) malloc(sizeof(particle_buffer_t));
                memset(g->outgoing_local[i].buffer, 0, sizeof(particle_buffer_t));
                g->outgoing_local[i].octant_id = new_id;
                g->outgoing_local[i].buffer->capacity = 10;
                 //10 particle space per local send buffer initally
                g->outgoing_local[i].buffer->particles = (particle_t *) malloc(g->outgoing_local[i].buffer->capacity * sizeof(particle_t));
                printf("\nTRIPLE sanity check local send");
                insert_particle(p, g->outgoing_local[i].buffer);
            }
            break;
        }
    } //insert into local using new_id
}

void remove_particle(particle_buffer_t * buffer, int index){
    if(buffer == NULL || index < 0 || index >= buffer->count){
        return;
    }

    if(index < buffer->count){
        buffer->particles[index] = buffer->particles[buffer->count - 1];
        memset(&buffer->particles[buffer->count - 1],/*set byte*/ 0, /*num bytes*/sizeof(particle_t));
    }
    buffer->count--;
}

void combine_local(local_data_t * g){
    if(g->outgoing_local == NULL) return;
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        //find corresponding incoming local buffer 
        for(int j = 0; j < g->mesh->local_num_quadrants; ++j){
            if(g->outgoing_local[j].octant_id == i){
                for(int k = 0; k < g->outgoing_local[j].buffer->count; ++k){
                    printf("\noct: %d, particle: %d being pushed-back", i, k);
                    insert_particle(&g->outgoing_local[j].buffer->particles[k], data->buffer);
                }
            }
        }
    }
}

void populate_send_buffers(local_data_t * g){
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
            p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
            octant_data_t * data = (octant_data_t *) oct->p.user_data;
            print_DEBUG_particle_data(g, data);
            for (int j = 0; j < data->buffer->count; ++j) {
                particle_t * p = &data->buffer->particles[j];

                int new_id = query_oct_id(g, p);

                //printf("oct: %d, paricle: %d, new_id:, %d", i, j, new_id);
                if(new_id == data->octant_id || /*DEBUG*/ new_id == -1) continue; //particle not transferred
                p8est_quadrant_t *oct_rec = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, new_id, NULL, NULL);
                octant_data_t * data2 = (octant_data_t *) oct_rec->p.user_data;
                //printf("\nbounds- x: %.2f, y:%.2f, z:%.2f \n", data2->bounds->x_min,data2->bounds->y_min, data2->bounds->z_min);
                insert_particle_into_outgoing(g, &data->buffer->particles[j], new_id);
                remove_particle(data->buffer, j--);
            }
        }
}

void do_dynamics(local_data_t * g){
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
            p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
            octant_data_t * data = (octant_data_t *) oct->p.user_data;
        for(int i = 0; i < data->buffer->count; i++){
            particle_t * p = &data->buffer->particles[i];

            //distance (negatives do matter)
            double dx = p->x - g->planet_xyz[0];
            double dy = p->y - g->planet_xyz[1];
            double dz = p->z - g->planet_xyz[2];
            double r_squared = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r_squared);
            double r3 = r_squared * r;

            //no div/0
            if (r3 < 1e-12) continue;

            double acc = -g->grav_const * g->planet_mass / r3;
            double ax = acc * dx;
            double ay = acc * dy;
            double az = acc * dz;

            //velo update
            p->vx += ax * g->timestep;
            p->vy += ay * g->timestep;
            p->vz += az * g->timestep;

            //pos update
            p->x += p->vx * g->timestep;
            p->y += p->vy * g->timestep;
            p->z += p->vz * g->timestep;
        }
    }
}

void loop(local_data_t * g){
    for(int cur_step = 0; cur_step < g->num_steps; ++cur_step){
        print_DEBUG_send_data(g);


        do_dynamics(g);
        
        populate_send_buffers(g);

        //Outgoing data has been prepared

        //combine local transfer data
        combine_local(g);

        printf("GOT TO SECOND ITERATION");
    }
}

void run(int argc, char **argv){

    local_data_t global, *g = &global; //rank-local data
    memset (g, 0, sizeof (*g));

    g->ppq = 2; //particles per octant
    g->planet_xyz[0] = 1.5; // mid = (1.5)^3 for a (3.0)^3 box
    g->planet_xyz[1] = 1.5;
    g->planet_xyz[2] = 1.5;
    g->planet_mass = 0.07;
    g->grav_const = 1.0;
    g->num_octants = 0;
    g->num_steps = 10;
    g->timestep = 1.0;

    mpi_context_t mpi_context;
    g->mpi = &mpi_context;
    mpi_set(argc, argv, g);

    p8est_setup(g);
    p8est_partition(g->p8est, 0, NULL); // 0 and NULL for uniform weight distribution

    g->ghost = p8est_ghost_new(g->p8est, P8EST_CONNECT_FULL);
    g->mesh = p8est_mesh_new(g->p8est, g->ghost, P8EST_CONNECT_FULL);

    loop(g);

    cleanup(g);
}

void print_DEBUG_particle_data(const local_data_t * g, const octant_data_t * data){
    printf("\n\n");
	for(int i = 0; i < data->buffer->count; i++){
		if(i == 0){
		       	printf("\n\n Particle data of octant:%d -- \n", data->octant_id);
			printf("Octant Bounds: x: %.1f-%.1f, y: %.1f-%.1f, z:%.1f-%.1f",
				       	data->bounds->x_min,
				       	data->bounds->x_max, 
					data->bounds->y_min,
				       	data->bounds->y_max, 
					data->bounds->z_min, 
					data->bounds->z_max);
		}
		printf("\n\nparticle :%d x: %.3f, y: %.3f, z: %.3f", i, data->buffer->particles[i].x, data->buffer->particles[i].y, data->buffer->particles[i].z);
        printf("\n\nparticle :%d vx: %.3f, vy: %.3f, vz: %.3f", i, data->buffer->particles[i].vx, data->buffer->particles[i].vy, data->buffer->particles[i].vz);

		if(i == data->buffer->count - 1) printf("\n\n");
	}
}
void print_DEBUG_send_data(const local_data_t * g){
    if(g->outgoing_local == NULL) printf("NO LOCAL SEND");
    //if(outgoing_ghost == NULL) printf("NO GHOST SEND");

    if(g->outgoing_local != NULL){
        for(int i = 0; i < g->mesh->local_num_quadrants; ++i){
            octant_data_t * data = &g->outgoing_local[i];
            printf("\n\nDATA going to octant %d\n", i);
            if(data->buffer && data->buffer->particles){
                for(int j = 0; j < data->buffer->count; ++j){
                    printf("particle %d, x: %.3f, y:%.3f, z:%.3f\n", j, data->buffer->particles[j].x, data->buffer->particles[j].y, data->buffer->particles[j].z);
                }
            }
        }
    }
}

int main(int argc, char **argv){
    run(argc, argv);
    return 0;
}
