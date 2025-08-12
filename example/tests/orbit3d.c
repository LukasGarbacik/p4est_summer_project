#include "orbit3d.h"

static sc_rand_state_t rand_object;
static int num_switch_id = 0;
static int num_insert = 0;
void free_particle_data(local_data_t * g){
    if(!g || !g->mesh) return;
    
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        if(!oct) continue;
        
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        if(data){
            if(data->buffer){
                if(data->buffer->particles) {
                    free(data->buffer->particles);
                    data->buffer->particles = NULL;
                }
                free(data->buffer);
                data->buffer = NULL;
            }
            free(data);
            oct->p.user_data = NULL;
        }
    }
}

void cleanup(local_data_t * g){
    if(!g) return;
    if(g->ghost_data) {
        sc_array_destroy(g->ghost_data);
        free(g->ghost_data);
        g->ghost_data = NULL;
    }
    
    //free_particle_data(g);
    //DEBUG
    //free count **'s
    
    if(g->mesh) {
        p8est_mesh_destroy(g->mesh);
        g->mesh = NULL;
    }
    if(g->ghost) {
        p8est_ghost_destroy(g->ghost);
        g->ghost = NULL;
    }
    if(g->p8est) {
        p8est_destroy(g->p8est);
        g->p8est = NULL;
    }
    if(g->connectivity) {
        p8est_connectivity_destroy(g->connectivity);
        g->connectivity = NULL;
    }
    
    if(g->mpi) {
        sc_MPI_Barrier(g->mpi->mpicomm);
        int mpiret = sc_MPI_Finalize();
        SC_CHECK_MPI(mpiret);
    }
}

void mpi_set(int argc, char **argv, local_data_t * data) {
    if(!data || !data->mpi) return;
    
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);

    data->mpi->mpicomm = sc_MPI_COMM_WORLD;
    sc_MPI_Comm_rank(data->mpi->mpicomm, &data->mpi->mpirank);
    sc_MPI_Comm_size(data->mpi->mpicomm, &data->mpi->mpisize);

    sc_init(data->mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
}

void p8est_setup(local_data_t * data){
    if(!data) return;
    
    p4est_init(NULL, SC_LP_DEFAULT);

    data->connectivity = p8est_connectivity_new_brick(3, 3, 3, 1, 1, 1);
    if(!data->connectivity) {
        printf("ERROR: Failed to create connectivity\n");
        return;
    }

    data->p8est = p8est_new_ext(data->mpi->mpicomm, data->connectivity, 
                        0,//Min number of quadrants
                        0,//refine level
                        1,//uniform t/f
                        sizeof(octant_data_t),//size of octant data struct
                        oct_init,//quadrant initialization function
                        data);//pointer for global data struct
    
    if(!data->p8est) {
        printf("ERROR: Failed to create p8est\n");
        return;
    }
    
    data->p8est->user_pointer = data;
}

void populate_oct_bounds(local_data_t * g, p4est_topidx_t which_tree, p8est_quadrant_t *octant){
    if(!g || !octant || !octant->p.user_data) return;
    
    octant_data_t *data = (octant_data_t *) octant->p.user_data;

    int x = octant->x;
    int y = octant->y;
    int z = octant->z;

    double vertex[3];
    p8est_qcoord_to_vertex(g->connectivity, which_tree, x, y, z, vertex);

    data->bounds.x_min = vertex[0];
    data->bounds.x_max = vertex[0] + 1.0;
    data->bounds.y_min = vertex[1];
    data->bounds.y_max = vertex[1] + 1.0;
    data->bounds.z_min = vertex[2];
    data->bounds.z_max = vertex[2] + 1.0;
}

void populate_particles(local_data_t * g, octant_data_t * data){
    if(!g || !data || !data->buffer || !data->buffer->particles) return;
    
    // Initialize random state with rank to ensure different seeds per process
    rand_object = (uint64_t)(time(NULL) + g->mpi->mpirank * 1000);
    
    for(int i = 0; i < data->buffer->count; i++){
        particle_t * p = &data->buffer->particles[i];

        p->octant_id = data->octant_id;
        
        // Generate random positions within bounds
        double rand_x = sc_rand(&rand_object);
        double rand_y = sc_rand(&rand_object);
        double rand_z = sc_rand(&rand_object);
        
        p->x = data->bounds.x_min + rand_x * (data->bounds.x_max - data->bounds.x_min); //always 1 for max-min but could be scaled different for other connectivities
        p->y = data->bounds.y_min + rand_y * (data->bounds.y_max - data->bounds.y_min);
        p->z = data->bounds.z_min + rand_z * (data->bounds.z_max - data->bounds.z_min);

        //calculate tangential velocity in 3D
        double dx = p->x - g->planet_xyz[0];
        double dy = p->y - g->planet_xyz[1];
        double dz = p->z - g->planet_xyz[2];
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        
        if(r < 1e-12) {
            printf("WARNING: Particle too close to planet center\n");
            p->vx = p->vy = p->vz = 0.0;
            continue;
        }
        
        double v = sqrt(g->grav_const * g->planet_mass / r);
        
        //Arbitrary vector
        double ax = 0.7, ay = 1.0, az = 0.4;
        if (fabs(ax - dx) < 1e-6 && fabs(ay - dy) < 1e-6) {
            ax = 0.2; ay = 0.5; az = 1.0;
        }
        
        //Xross product
        double x_vec = dy * az - dz * ay;
        double y_vec = dz * ax - dx * az;
        double z_vec = dx * ay - dy * ax;
        double mag = sqrt(x_vec*x_vec + y_vec*y_vec + z_vec*z_vec);
        
        if(mag < 1e-12) {
            printf("WARNING: Zero tangential velocity\n");
            p->vx = p->vy = p->vz = 0.0;
            continue;
        }

        p->vx = v * x_vec / mag;
        p->vy = v * y_vec / mag;
        p->vz = v * z_vec / mag;

        if(!in_bounds(p, &data->bounds)) {
            //printf("ERROR: Particle initialized outside bounds\n");
        }
    }
}

void oct_init(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *octant){
    if(!p8est || !octant) return;
    
    local_data_t *g = (local_data_t *) p8est->user_pointer;
    if(!g) return;
    
    octant->p.user_data = calloc(1, sizeof(octant_data_t));
    if(!octant->p.user_data) {
        printf("ERROR: Failed to allocate octant data\n");
        return;
    }
    
    octant_data_t *data = (octant_data_t *) octant->p.user_data;
    data->buffer = (particle_buffer_t *) calloc(1, sizeof(particle_buffer_t));
    if(!data->buffer) {
        printf("ERROR: Failed to allocate particle buffer\n");
        free(octant->p.user_data);
        octant->p.user_data = NULL;
        return;
    }
    data->buffer->count = g->ppq;
    data->buffer->capacity = g->ppq * 2;
    data->buffer->particles = (particle_t *) calloc(data->buffer->capacity, sizeof(particle_t));
    if(!data->buffer->particles) {
        printf("ERROR: Failed to allocate particles array\n");
        free(data->buffer);
        free(octant->p.user_data);
        octant->p.user_data = NULL;
        return;
    }
    data->octant_id = g->num_octants++;
    data->mpirank = g->mpi->mpirank;

    populate_oct_bounds(g, which_tree, octant);
    populate_particles(g, data);
}

bool in_bounds(particle_t * p, octant_bounds_t * bounds){
    if(!p || !bounds) return false;
    
    return (p->x >= bounds->x_min && p->x <= bounds->x_max &&
            p->y >= bounds->y_min && p->y <= bounds->y_max &&
            p->z >= bounds->z_min && p->z <= bounds->z_max);
}

int query_oct_id(local_data_t * g, particle_t * p){
    if(!g || !p || !g->mesh) return -1;
    
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        if(!oct || !oct->p.user_data) continue;
        
        octant_data_t * data = (octant_data_t *) oct->p.user_data;

        if(in_bounds(p, &data->bounds)){
            //if(g->mpi->mpirank == 0) //printf("\nnew oct: %d, particle (%.2f, %.2f, %.2f) being sent", data->octant_id, p->x, p->y, p->z);
            //if(g->mpi->mpirank == 0) //printf("\nbounds- x: %.2f, y:%.2f, z:%.2f \n", data->bounds.x_min, data->bounds.y_min, data->bounds.z_min);
            return data->octant_id;
        }
    }
    return -1;//ghost
}

void insert_particle(particle_t * p, particle_buffer_t * buffer){
    if(!p || !buffer) {
        printf("ERROR: Null pointer in insert_particle\n");
        return;
    }
    
    num_insert++;
    
    if(!buffer->particles){
        printf("ERROR: Null particles array in insert_particle (INSERT NUMBER %d)\n", num_insert);
        return;
    }
    
    //dynamic reallocation
    if(buffer->count >= buffer->capacity){
        size_t new_capacity = buffer->capacity * 2;
        particle_t* new_particles = realloc(buffer->particles, new_capacity * sizeof(particle_t));
        if(!new_particles){
            printf("ERROR: Failed to reallocate particles buffer from %zu to %zu\n", 
                   buffer->capacity, new_capacity);
            return;
        }
        buffer->particles = new_particles;
        buffer->capacity = new_capacity;
    }
    
    
    buffer->particles[buffer->count] = *p;
    buffer->count++;
}

void insert_particle_into_outgoing(local_data_t * g, p8est_quadrant_t *oct, particle_t * p, int new_id){
    if(!g || !p) return;
    
    if(new_id == -1){
        insert_ghost_particle(g, oct, p);
        return;
    }
    
    
    if(!g->outgoing_local){
        g->outgoing_local = (octant_data_t *) calloc(g->mesh->local_num_quadrants, sizeof(octant_data_t));
        if(!g->outgoing_local) {
            printf("ERROR: Failed to allocate outgoing_local\n");
            return;
        }
        
        for(int i = 0; i < g->mesh->local_num_quadrants; ++i){
            g->outgoing_local[i].octant_id = -2; // dummy value
            g->outgoing_local[i].buffer = NULL;
        }
    }
    
    // Find existing buffer or create new one
    for(int i = 0; i < g->mesh->local_num_quadrants; ++i){
        if(g->outgoing_local[i].octant_id == new_id){
            // Found existing buffer
            insert_particle(p, g->outgoing_local[i].buffer);
            return;
        }
        else if(g->outgoing_local[i].octant_id == -2){
            // First unused slot - initialize it
            g->outgoing_local[i].buffer = (particle_buffer_t *) calloc(1, sizeof(particle_buffer_t));
            if(!g->outgoing_local[i].buffer) {
                printf("ERROR: Failed to allocate particle buffer\n");
                return;
            }
            
            g->outgoing_local[i].octant_id = new_id;
            g->outgoing_local[i].buffer->capacity = 10;
            g->outgoing_local[i].buffer->count = 0;
            g->outgoing_local[i].buffer->particles = (particle_t *) calloc(10, sizeof(particle_t));
            
            if(!g->outgoing_local[i].buffer->particles) {
                printf("ERROR: Failed to allocate particles array\n");
                free(g->outgoing_local[i].buffer);
                g->outgoing_local[i].buffer = NULL;
                return;
            }
            
            insert_particle(p, g->outgoing_local[i].buffer);
            return;
        }
    }
    printf("ERROR: No available slot for octant %d\n", new_id);
}

void remove_particle(particle_buffer_t * buffer, int index){
    if(!buffer || !buffer->particles || index < 0 || index >= buffer->count){
        printf("ERROR: Invalid remove_particle call - buffer=%p, index=%d, count=%d\n", 
               (void*)buffer, index, buffer ? buffer->count : -1);
        return;
    }

    if(index < buffer->count - 1){
        buffer->particles[index] = buffer->particles[buffer->count - 1];
    }
    
    memset(&buffer->particles[buffer->count - 1], 0, sizeof(particle_t)); //if removing the last, then this clears it
    buffer->count--;
}

void free_local(local_data_t * g){
    if(!g || !g->outgoing_local) return;
    
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        if(g->outgoing_local[i].buffer){
            if(g->outgoing_local[i].buffer->particles) {
                free(g->outgoing_local[i].buffer->particles);
                g->outgoing_local[i].buffer->particles = NULL;
            }
            free(g->outgoing_local[i].buffer);
            g->outgoing_local[i].buffer = NULL;
        }
    }
    free(g->outgoing_local);
    g->outgoing_local = NULL;
}

void combine_local(local_data_t * g){
    if(!g || !g->outgoing_local || !g->mesh) return;
    
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        if(!oct || !oct->p.user_data) continue;
        
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        if(!data->buffer) continue;
        
        //find corresponding incoming local buffer 
        for(int j = 0; j < g->mesh->local_num_quadrants; ++j){
            if(g->outgoing_local[j].octant_id == data->octant_id && 
               g->outgoing_local[j].buffer && 
               g->outgoing_local[j].buffer->particles){
                
                for(int k = 0; k < g->outgoing_local[j].buffer->count; ++k){
                    insert_particle(&g->outgoing_local[j].buffer->particles[k], data->buffer);
                }
                break;
            }
        }
    }
}

void populate_send_buffers(local_data_t * g){
    if(!g || !g->mesh) return;
    if(g->mpi->mpisize > 1 && g->outgoing_ghost == NULL){//always initalized even before insertion (allow recv ghost on any loop)
        g->outgoing_ghost = (octant_data_t **) calloc(g->mpi->mpisize, sizeof(octant_data_t *));
        if(!g->outgoing_ghost){
            printf("\nFailed to allocate ghost 2d array");
            return;
        }
    }
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        if(!oct || !oct->p.user_data) continue;
        
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        if(!data->buffer || !data->buffer->particles) continue;
        
        
        //Process backwards
        for (int j = data->buffer->count - 1; j >= 0; --j) {
            particle_t * p = &data->buffer->particles[j];
            int new_id = query_oct_id(g, p);

            if(new_id == data->octant_id) continue; // particle not transferred
            
            num_switch_id++;
            insert_particle_into_outgoing(g, oct, &data->buffer->particles[j], new_id);
            remove_particle(data->buffer, j);
        }
    }
}

void do_dynamics(local_data_t * g){
    if(!g || !g->mesh) return;
    
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        if(!oct || !oct->p.user_data) continue;
        octant_data_t * data = (octant_data_t *) oct->p.user_data;

        
        if(!data->buffer || !data->buffer->particles) continue;
        
        for(int j = 0; j < data->buffer->count; j++){
            particle_t * p = &data->buffer->particles[j];

            double dx = p->x - g->planet_xyz[0];
            double dy = p->y - g->planet_xyz[1];
            double dz = p->z - g->planet_xyz[2];
            double r_squared = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r_squared);
            double r3 = r_squared * r;

            //division by zero
            if (r3 < 1e-12) continue;

            double acc = -g->grav_const * g->planet_mass / r3;
            double ax = acc * dx;
            double ay = acc * dy;
            double az = acc * dz;

            p->vx += ax * g->timestep;
            p->vy += ay * g->timestep;
            p->vz += az * g->timestep;

            p->x += p->vx * g->timestep;
            p->y += p->vy * g->timestep;
            p->z += p->vz * g->timestep;

            //global wrap
            double sum_trick = p->x + p->y + p->z; //less lines than keeping a bool checker
            if(p->x > global_bound) p->x -= global_bound;
            if(p->x < 0) p->x += global_bound;

            if(p->y > global_bound) p->y -= global_bound;
            if(p->y < 0) p->y += global_bound;

            if(p->z > global_bound) p->z -= global_bound;
            if(p->z < 0) p->z += global_bound;

            if(p->x + p->y + p->z != sum_trick){//retangitize after periodic switch
                double dx = p->x - g->planet_xyz[0];
                double dy = p->y - g->planet_xyz[1];
                double dz = p->z - g->planet_xyz[2];
                double r = sqrt(dx*dx + dy*dy + dz*dz);

                double v = sqrt(g->grav_const * g->planet_mass / r);
        
                //Arbitrary vector
                double ax = 0.7, ay = 1.0, az = 0.4;
                if (fabs(ax - dx) < 1e-6 && fabs(ay - dy) < 1e-6) {
                    ax = 0.2; ay = 0.5; az = 1.0;
                }

                //Xross product
                double x_vec = dy * az - dz * ay;
                double y_vec = dz * ax - dx * az;
                double z_vec = dx * ay - dy * ax;
                double mag = sqrt(x_vec*x_vec + y_vec*y_vec + z_vec*z_vec);

                if(mag < 1e-12) {
                    //printf("WARNING: Zero tangential velocity\n");
                    p->vx = p->vy = p->vz = 0.0;
                    continue;
                }

                p->vx = v * x_vec / mag;
                p->vy = v * y_vec / mag;
                p->vz = v * z_vec / mag;
            }
        }
    }
}

void insert_ghost_particle(local_data_t * g, p8est_quadrant_t *oct, particle_t * p){
    octant_data_t * data = (octant_data_t *) oct->p.user_data;
    //g and p are always valid from insert_particle_into_outgoing
    for (int i = 0; i < g->ghost->ghosts.elem_count; ++i) {
        octant_data_t * ghost_data = (octant_data_t *) sc_array_index(g->ghost_data, i);
        if(!ghost_data){
            printf("ERROR: GHOST OCTANT DATA IS NOT ALLOCATED");
            return;
        }
        if(in_bounds(p, &ghost_data->bounds)){
            p->octant_id = ghost_data->octant_id;//new octant id for outgoing particle
            //prep to put particle p in g->outgoing_ghost[ghost_data->mpirank]
            if(!g->outgoing_ghost[ghost_data->mpirank]){
                g->outgoing_ghost[ghost_data->mpirank] = (octant_data_t *) calloc(g->no, sizeof(octant_data_t));
            }
            //index specific octant too
            octant_data_t * insert_oct = &g->outgoing_ghost[ghost_data->mpirank][ghost_data->octant_id];
            if(!insert_oct->buffer){
                //match octant id to ghost access so it can be easily combined after mpi
                insert_oct->octant_id = ghost_data->octant_id;
                insert_oct->buffer = (particle_buffer_t *) calloc(1, sizeof(particle_buffer_t));
                insert_oct->buffer->capacity = 10;
                insert_oct->buffer->count = 0;
                insert_oct->buffer->particles = (particle_t *) calloc(10, sizeof(particle_t));
            }
            insert_particle(p, insert_oct->buffer);
        }
    }
}

void free_ghost(local_data_t * g){
    if(!g->outgoing_ghost) return;
    for(int i = 0; i < g->mpi->mpisize; i++){
        if(!g->outgoing_ghost[i]) continue;
        for(int j = 0; j < g->no; j++){
            if(g->outgoing_ghost[i][j].buffer){
                if(g->outgoing_ghost[i][j].buffer->particles){
                    free(g->outgoing_ghost[i][j].buffer->particles);
                    g->outgoing_ghost[i][j].buffer->particles = NULL;
                }
                free(g->outgoing_ghost[i][j].buffer);
                g->outgoing_ghost[i][j].buffer = NULL;
            }
        }
        free(g->outgoing_ghost[i]);
        g->outgoing_ghost[i] = NULL;
    }

    for(int i = 0; i < g->mpi->mpisize; i++){
        if(g->send_raw_data && g->send_raw_data[i]){
            free(g->send_raw_data[i]);
            g->send_raw_data[i] = NULL;
        }
        if(g->recv_raw_data && g->recv_raw_data[i]){
            free(g->recv_raw_data[i]);
            g->recv_raw_data[i] = NULL;
        }
    }

    free(g->outgoing_ghost);
    g->outgoing_ghost = NULL;
    free(g->send_raw_data);
    g->send_raw_data = NULL;
    free(g->recv_raw_data);
    g->recv_raw_data = NULL;
    free(g->send_bytes_count);
    g->send_bytes_count = NULL;
    free(g->recv_bytes_count);
    g->recv_bytes_count = NULL;
}

void populate_byte_buffers(local_data_t * g){
    if(!g->outgoing_ghost) return;
    if(!g->send_bytes_count) g->send_bytes_count = (int *) calloc(g->mpi->mpisize, sizeof(int));
    if(!g->recv_bytes_count) g->recv_bytes_count = (int *) calloc(g->mpi->mpisize, sizeof(int));

    //if given rank or octant is not allocated, its array(s) will be 0 for bytes
    
    //send_bytes_count[i] format:
    //int #octants
    //loop:
    //local_octant_id
    //^count
    //particle raw data
    for(int rank = 0; rank < g->mpi->mpisize; rank++){
        if(!g->outgoing_ghost[rank]){
            g->send_bytes_count[rank] = 0;
            continue;   
        }
        int total = 0;
        total += sizeof(int); //#octants
        for(int oct = 0; oct < g->no; oct++){
            if(g->outgoing_ghost[rank][oct].buffer){
                total += 2 * sizeof(int); //local oct_id and particle count
                total += g->outgoing_ghost[rank][oct].buffer->count * sizeof(particle_t);
            }
        }
        g->send_bytes_count[rank] = total;
    }

    if(!g->send_raw_data){
        g->send_raw_data = (char **) calloc(g->mpi->mpisize, sizeof(char *));
    }
    if(!g->recv_raw_data){
        g->recv_raw_data = (char **) calloc(g->mpi->mpisize, sizeof(char *));
    }

    for(int rank = 0; rank < g->mpi->mpisize; rank++){
        g->send_raw_data[rank] = (char *) malloc(g->send_bytes_count[rank]); //can do malloc(0), and free, but dont access
        char * ptr = g->send_raw_data[rank];

        if(!g->outgoing_ghost[rank]) continue;

        int num_oct_count = 0;
        ptr += sizeof(int);//skip total oct until all loops finish

        for(int oct = 0; oct < g->no; oct++){
            if(g->outgoing_ghost[rank][oct].buffer){
                num_oct_count++;
                if(!g->outgoing_ghost[rank][oct].buffer->particles){
                    printf("\nERROR: particle buffer wrongfully not allocated when doing byte conversion");
                    return;
                }

                //copy octant id to buffer
                memcpy(ptr, &g->outgoing_ghost[rank][oct].octant_id, sizeof(int));
                ptr += sizeof(int);

                //copy given oct's buffer count id to buffer
                int count = g->outgoing_ghost[rank][oct].buffer->count;
                memcpy(ptr, &count, sizeof(int));
                ptr += sizeof(int);
                
                //copy given particle buffer into byte buffer
                memcpy(ptr, g->outgoing_ghost[rank][oct].buffer->particles, count * sizeof(particle_t));
                ptr += count * sizeof(particle_t);
            }
        }
        memcpy(g->send_raw_data[rank], &num_oct_count, sizeof(int));
    }
}


void send_recv_counts(local_data_t * g){
    if(g->outgoing_ghost == NULL)  return;


    int index = 0;//accounting for the self skip on rank == rank
    MPI_Request recieved_promises[g->mpi->mpisize - 1];
    MPI_Request sent_promises[g->mpi->mpisize - 1];

    //recv bytes
    for(int rank = 0; rank < g->mpi->mpisize; rank++){
        if(g->mpi->mpirank == rank) continue;
        MPI_Irecv(
                /*local memory*/&g->recv_bytes_count[rank],
                /* count */1,
                MPI_INT, 
                /* from rank:*/rank, 
                /* tag */ 0,
                g->mpi->mpicomm,
                &recieved_promises[index]);
        index++;
    }
    //send bytes
    index = 0;
    for(int rank = 0; rank < g->mpi->mpisize; ++rank){
        if(g->mpi->mpirank == rank) continue;
        MPI_Isend(
                /*local memory*/&g->send_bytes_count[rank],
                /* count */1,
                MPI_INT, 
                /*  rank*/rank, 
                /* tag */ 0,
                g->mpi->mpicomm,
                &sent_promises[index]);
        index++;
    }
    MPI_Waitall(g->mpi->mpisize - 1, recieved_promises, MPI_STATUSES_IGNORE);
    MPI_Waitall(g->mpi->mpisize - 1, sent_promises, MPI_STATUSES_IGNORE);
}

void send_recv_particles(local_data_t * g){
    if(g->send_bytes_count == NULL || g->recv_bytes_count == NULL) return;
    if(g->send_raw_data == NULL || g->recv_raw_data == NULL) return;
    if(g->outgoing_ghost == NULL)  return;


    //allocate temp data for direct mpi recv based on the buffer sizes in recv_bytes_count
    for(int rank = 0; rank < g->mpi->mpisize; rank++){
        if(g->recv_bytes_count[rank] > 0){
            g->recv_raw_data[rank] = (char *) malloc(g->recv_bytes_count[rank]);
        }
    }


    int index = 0;//accounting for the self skip on rank == rank
    MPI_Request recieved_promises[g->mpi->mpisize - 1];
    MPI_Request sent_promises[g->mpi->mpisize - 1];

    //if given send is 0, skip send and fill promise on both sides
    //else, send prepared buffer

    //recv raw buffers
    for(int rank = 0; rank < g->mpi->mpisize; rank++){
        if(g->mpi->mpirank == rank) continue;
        if(g->recv_bytes_count[rank] == 0){//dummy recv for 0 count
            MPI_Irecv(NULL, 0, MPI_BYTE, rank, 1, g->mpi->mpicomm, &recieved_promises[index]);
            index++;
            continue;
        }
        MPI_Irecv(
            /*local memory*/g->recv_raw_data[rank],
            /* count */g->recv_bytes_count[rank],
            MPI_BYTE,
            rank,
            1, //tag
            g->mpi->mpicomm,
            &recieved_promises[index]
        );
        index++;
    }

    //send raw buffers
    index = 0; //reset index for send loop
    for(int rank = 0; rank < g->mpi->mpisize; rank++){
        if(g->mpi->mpirank == rank) continue;
        if(g->send_bytes_count[rank] == 0){//dummy recv for 0 count
            MPI_Isend(NULL, 0, MPI_BYTE, rank, 1, g->mpi->mpicomm, &sent_promises[index]);
            index++;
            continue;
        }
        MPI_Isend(
                /*local memory*/g->send_raw_data[rank],
                /* count */g->send_bytes_count[rank],
                MPI_BYTE, 
                /*  rank*/rank, 
                /* tag */ 1,
                g->mpi->mpicomm,
                &sent_promises[index]
        );
        index++;
    }
    MPI_Waitall(g->mpi->mpisize - 1, recieved_promises, MPI_STATUSES_IGNORE);
    MPI_Waitall(g->mpi->mpisize - 1, sent_promises, MPI_STATUSES_IGNORE);
}

void combine_ghost(local_data_t * g){
    if(g->outgoing_ghost == NULL)  return;
    if(g->recv_raw_data == NULL || g->recv_bytes_count == NULL) return;

    //dynamic inner loop with resetting buffer pointer relative to the count of octs at the front
    for(int rank = 0; rank < g->mpi->mpisize; rank++){
        if(g->mpi->mpirank == rank || g->recv_raw_data[rank] == NULL) continue;
        char * ptr = g->recv_raw_data[rank];
        int incoming_octants;
        memcpy(&incoming_octants, ptr, sizeof(int));
        ptr += sizeof(int);
        for(int oct_index = 0; oct_index < incoming_octants; oct_index++){
            int insert_id;
            memcpy(&insert_id, ptr, sizeof(int));
            ptr += sizeof(int);
            if(insert_id >= g->mesh->local_num_quadrants){
                printf("\nERROR: Invalid incoming octant id");
                return;
            }
            p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, insert_id, NULL, NULL);
            if(!oct || !oct->p.user_data) continue;
            octant_data_t * data = (octant_data_t *) oct->p.user_data;
            if(data->buffer == NULL || data->buffer->particles == NULL){
                printf("\nLocal octant data wrongfully not allocated when combining ghost");
                return;
            }

            int incoming_count;
            memcpy(&incoming_count, ptr, sizeof(int));
            ptr += sizeof(int);

            //copy (incoming_count) # particles from ptr to data->buffer->particles
            for(int i = 0; i < incoming_count; i++){
                particle_t particle;
                memcpy(&particle, ptr, sizeof(particle_t));
                ptr += sizeof(particle_t);

                if(g->mpi->mpirank == 0){
                    printf("\nParticle info: id:%d x:%.3f, y:%.3f, z:%.3f", particle.octant_id, particle.x, particle.y, particle.z);
                }

                insert_particle(&particle, data->buffer);
            }
        }
    }
}

void write_vtk(local_data_t * g, const char *output_dir, int cur_step){
    char vtk_filename[256];
    snprintf(vtk_filename, sizeof(vtk_filename), "%s/particles_rank%d_step%d.vtk", output_dir, g->mpi->mpirank, cur_step);
    FILE *f = fopen(vtk_filename, "w");
    if (!f) {
        printf("Error opening %s\n", vtk_filename);
        return;
    }
    int total_particles = 0;
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        total_particles += data->buffer->count;
    }

    fprintf(f, "# vtk DataFile Version 3.0\nParticles\nASCII\nDATASET POLYDATA\n");
    fprintf(f, "POINTS %d float\n", total_particles);
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        for (int j = 0; j < data->buffer->count; ++j) {
            fprintf(f, "%f %f %f\n", data->buffer->particles[j].x, data->buffer->particles[j].y, data->buffer->particles[j].z);
        }
    }
    fprintf(f, "VERTICES %d %d\n", total_particles, 2 * total_particles);
    for (int i = 0; i < total_particles; ++i) {
        fprintf(f, "1 %d\n", i);
    }
    fprintf(f, "POINT_DATA %d\n", total_particles);
    fprintf(f, "VECTORS velocity float\n");
    for (p4est_locidx_t i = 0; i < g->mesh->local_num_quadrants; ++i) {
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        for (int j = 0; j < data->buffer->count; ++j) {
            fprintf(f, "%f %f %f\n", data->buffer->particles[j].vx, data->buffer->particles[j].vy, data->buffer->particles[j].vz);
        }
    }
    fclose(f);
}


void loop(local_data_t * g, const char *output_dir){
    if(!g) return;
    
    for(int cur_step = 0; cur_step < g->num_steps; ++cur_step){
        //visual output
        write_vtk(g, output_dir, cur_step);
        
        //single step dynamics
        do_dynamics(g);
        
        //g->outgoing_local/ghost
        populate_send_buffers(g);

        //debug output
        print_long_DEBUG(g, cur_step);
        
        //combine local data
        combine_local(g);
        
        //prep ghost mpi transfer
        populate_byte_buffers(g);
        send_recv_counts(g);
        
        //ghost transfer
        send_recv_particles(g);

        //merge incoming with local
        combine_ghost(g);

        //clear memory
        free_local(g);
        free_ghost(g);

        int total = 0;
        for(int i = 0; i < g->mesh->local_num_quadrants; i++){
            p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
            if(!oct || !oct->p.user_data) continue;
            octant_data_t * data = (octant_data_t *) oct->p.user_data;
            total += data->buffer->count;
        }
        //printf("\nRANK:%d has a total of %d particles", g->mpi->mpirank, total);
    }
}
void print_long_DEBUG(local_data_t * g, int loop_num){
    char vtk_filename[256];
    snprintf(vtk_filename, sizeof(vtk_filename), "debug_output_rank%d.txt", g->mpi->mpirank);
    FILE *f = fopen(vtk_filename, "a");
    if (!f) {
        printf("Error opening %s\n", vtk_filename);
        return;
    }

    fprintf(f,"\n\nLoop number:%d", loop_num);
    fflush(f);
    fprintf(f,"\nLocal particle data");
    fflush(f);
    for(int i = 0; i < g->mesh->local_num_quadrants; i++){
        p8est_quadrant_t *oct = p8est_mesh_quadrant_cumulative(g->p8est, g->mesh, i, NULL, NULL);
        octant_data_t * data = (octant_data_t *) oct->p.user_data;
        if(!data){
            fprintf(f,"\nNO OCT DATA");
            fflush(f);
            return;
        }
        fprintf(f,"\n\nParticles for octant:%d", i);
        fflush(f);
        for(int j = 0; j < data->buffer->count; j++){
            fprintf(f,"\nparticle %d: x:%.3f, y:%.3f, z:%.3f", j, data->buffer->particles[j].x, data->buffer->particles[j].y, data->buffer->particles[j].z);
            fflush(f);
        }
    }
    fprintf(f,"\n\nOutgoing local data\n");
    fflush(f);
    if(!g->outgoing_local){
        fprintf(f,"\nNo outgoing local data this loop");
        fflush(f);
    }
    else{
        for(int i = 0; i < g->mesh->local_num_quadrants; i++){
            octant_data_t * data = &g->outgoing_local[i];
            if(!data) continue;
            if(data->octant_id < 0) continue;
            fprintf(f,"\ndata going to octant: %d", data->octant_id);
            fflush(f);
            for(int j = 0; j < data->buffer->count; j++){
                fprintf(f,"\nparticle %d: x:%.3f, y:%.3f, z:%.3f", j, data->buffer->particles[j].x, data->buffer->particles[j].y, data->buffer->particles[j].z);
                fflush(f);
            }
        }
    }
    //outgoing ghost buffers
    fprintf(f,"\n\nOutgoing ghost data\n");
    fflush(f);
    bool ghost = false;
    if(g->outgoing_ghost){
        for(int i = 0; i < g->mpi->mpisize; i++){
            octant_data_t * rank_array = g->outgoing_ghost[i];
            if(!rank_array) continue;
            ghost = true;
            fprintf(f,"\ndata going to rank: %d", i);
            fflush(f);
            for(int j = 0; j < g->no; j++){
                if(!rank_array[j].buffer) continue;
                fprintf(f,"\ndata going to rank:%d, octant:%d", i, j);
                fflush(f);
                for(int k = 0; k < rank_array[j].buffer->count; k++){
                    fprintf(f,"\nparticle %d: x:%.3f, y:%.3f, z:%.3f", rank_array[j].buffer->particles[k].x, rank_array[j].buffer->particles[k].y, rank_array[j].buffer->particles[k].z);
                    fflush(f);
                }
            }
        }
    }
    if(!ghost){
        fprintf(f, "\nNo ghost data this loop");
        fflush(f);
    }
    fclose(f);
}

void run(int argc, char **argv){
    local_data_t global;
    local_data_t *g = &global;
    memset(g, 0, sizeof(*g));

    // Initialize parameters
    g->ppq = 7; // particles per octant
    g->planet_xyz[0] = 1.5; // mid = (1.5)^3 for a (3.0)^3 box
    g->planet_xyz[1] = 1.5;
    g->planet_xyz[2] = 1.5;
    g->planet_mass = 0.07;
    g->grav_const = 1.0;
    g->num_octants = 0;
    g->num_steps = 200;
    g->timestep = 0.5;

    mpi_context_t mpi_context;
    memset(&mpi_context, 0, sizeof(mpi_context));
    g->mpi = &mpi_context;
    
    mpi_set(argc, argv, g);
    p8est_setup(g);
    
    if(!g->p8est) {
        printf("ERROR: Failed to setup p8est\n");
        cleanup(g);
        return;
    }
    if(g->mpi->mpisize > total_octants){
        printf("Too many processors for the number of octants in this mesh (%d) > (%d)", g->mpi->mpisize, total_octants);
        return;
    }
    
    p8est_partition(g->p8est, 0, NULL); // 0 and NULL for uniform weight distribution

    g->ghost = p8est_ghost_new(g->p8est, P8EST_CONNECT_FULL);
    if(!g->ghost) {
        printf("ERROR: Failed to create ghost\n");
        cleanup(g);
        return;
    }
    g->mesh = p8est_mesh_new(g->p8est, g->ghost, P8EST_CONNECT_FULL);
    if(!g->mesh) {
        printf("ERROR: Failed to create mesh\n");
        cleanup(g);
        return;
    }

    //dependant helper data filled after mesh and ghost
    g->total_considered_oct = g->mesh->local_num_quadrants + g->ghost->ghosts.elem_count;
    g->no = (g->total_considered_oct / g->mpi->mpisize) + 3;

    char output_dir[256];
    snprintf(output_dir, sizeof(output_dir), "output_vtk/particle_output_rank%d", g->mpi->mpirank);
    if(g->mpi->mpirank == 0) {
        mkdir("output_vtk", 0777); // Only rank 0 creates the main directory
    }
    sc_MPI_Barrier(g->mpi->mpicomm);//wait for rank one before other ranks add dirs
    mkdir(output_dir, 0777);//full permissions

    g->ghost_data = (sc_array_t*) malloc(sizeof(sc_array_t));
    if(!g->ghost_data) {
        printf("ERROR: Failed to allocate ghost_data\n");
        cleanup(g);
        return;
    }
    sc_array_init(g->ghost_data, sizeof(octant_data_t));
    sc_array_resize(g->ghost_data, g->ghost->ghosts.elem_count);
    p8est_ghost_exchange_data(g->p8est, g->ghost, g->ghost_data->array);

    loop(g, output_dir);
    printf("\nCLEANUP FUNCTION ON RANK: %d\n", g->mpi->mpirank);
    cleanup(g);
}
/*
void print_DEBUG_particle_data(const octant_data_t * data){
    if(!data || !data->buffer || !data->buffer->particles) return;
    
    //printf("\n\n");
    for(int i = 0; i < data->buffer->count; i++){
        //if(i == 0){
            //printf("\n\n Particle data of octant:%d -- \n", data->octant_id);
            //printf("Octant Bounds: x: %.1f-%.1f, y: %.1f-%.1f, z:%.1f-%.1f",
                   //data->bounds.x_min, data->bounds.x_max, 
                   //data->bounds.y_min, data->bounds.y_max, 
                   //data->bounds.z_min, data->bounds.z_max);
        //}
        //printf("\n\nparticle :%d x: %.3f, y: %.3f, z: %.3f", 
               //i, data->buffer->particles[i].x, data->buffer->particles[i].y, data->buffer->particles[i].z);
        //printf("\n\nparticle :%d vx: %.3f, vy: %.3f, vz: %.3f", 
               //i, data->buffer->particles[i].vx, data->buffer->particles[i].vy, data->buffer->particles[i].vz);

        //if(i == data->buffer->count - 1) //printf("\n\n");
    }
}

void print_DEBUG_send_data(const local_data_t * g){
    if(!g) return;
    
    if(!g->outgoing_local) {
        //printf("NO LOCAL SEND\n");
        return;
    }

    for(int i = 0; i < g->mesh->local_num_quadrants; ++i){
        if(g->outgoing_local[i].octant_id == -2) continue; // Skip unused slots
        
        if(g->outgoing_local[i].buffer && g->outgoing_local[i].buffer->particles){
            for(int j = 0; j < g->outgoing_local[i].buffer->count; ++j){
                //if(g->mpi->mpirank == 0) //printf("particle %d, x: %.3f, y:%.3f, z:%.3f\n", 
                //j,
                //g->outgoing_local[i].buffer->particles[j].x,
                //g->outgoing_local[i].buffer->particles[j].y,
                //g->outgoing_local[i].buffer->particles[j].z);
            }
        }
    }
}
*/
int main(int argc, char **argv){
    run(argc, argv);
    return 0;
}