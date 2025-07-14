#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <sc_notify.h>
#include <sc_options.h>
#include <sc_random.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "orbit.h"

#define p_per_quad 1

#define grav_const 1

#define timestep 0.1

static sc_rand_state_t global_rand_state;
static int global_mpi_rank = -1;

static const double planet_xyz[2] = {.50, .50}; //one planet right at the center (0.5, 0.5)
static const double planet_mass = 0.0167;

mpi_context_t mpi_init(int argc, char **argv) {

    mpi_context_t mpi_context;

    int mpiret = sc_MPI_Init(&argc, &argv);

    mpi_context.mpicomm = sc_MPI_COMM_WORLD;
    sc_MPI_Comm_rank(mpi_context.mpicomm, &mpi_context.mpirank);
    sc_MPI_Comm_size(mpi_context.mpicomm, &mpi_context.mpisize);

    SC_CHECK_MPI (mpiret);

    sc_init (mpi_context.mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

    global_rand_state = (uint64_t)(time(NULL) + mpi_context.mpirank);
    global_mpi_rank = mpi_context.mpirank;

    return mpi_context;
}

p4est_t * p4est_setup(mpi_context_t *mpi_context) {
    p4est_init(NULL, SC_LP_DEFAULT);

    p4est_connectivity_t *connectivity = p4est_connectivity_new_unitsquare();

    p4est_t *p4est = p4est_new_ext(mpi_context->mpicomm, connectivity, 0, 1, 0,
                                    /*octant struct size*/ sizeof(particle_buffer_t), 
                                    /*quadrant init fxn */ quad_init,
                                    /*pointer for global data struct */ NULL);

    return p4est;
}

particle_t particle_single_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant){
    
    double xs, ys;

    if(quadrant->x == 0 && quadrant->y == 0){
        xs = 0;
        ys = 0;
    }
    else if(quadrant->x == 0){
        xs = 0;
        ys = 0.5;
    }
    else if(quadrant->y == 0){
        xs = 0.5;
        ys = 0;
    }
    else{
        xs = 0.5;
        ys = 0.5;
    }

    
    particle_t particle;

    particle.quadrant_id = global_mpi_rank;
    
    //random positions rand_x/y in [0, 0.5)
    double rand_x = sc_rand(&global_rand_state) / 2;
    double rand_y = sc_rand(&global_rand_state) / 2;
    
    particle.x = xs + rand_x;
    particle.y = ys + rand_y;
    
    //velocity of particle calculation
    double dx, dy;
    dx = particle.x - planet_xyz[0];
    dy = particle.y - planet_xyz[1];
    double r = sqrt(dx*dx + dy*dy); //dist from planet
    double v = sqrt(grav_const * planet_mass / r); //magnitude of velocity

    particle.vx = -v * dy / r;
    particle.vy = v * dx / r;
    
    return particle;
}

void quad_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant){
    particle_buffer_t *buf = (particle_buffer_t *) quadrant->p.user_data;
    buf->count = p_per_quad;
    buf->capacity = p_per_quad;
    buf->particles = (particle_t *) malloc(4 * buf->count * sizeof(particle_t)); //4x for all particles (not dynamic)
    for(int i = 0; i < p_per_quad; i++){
        buf->particles[i] = particle_single_init(p4est, which_tree, quadrant);
    }
}

int find_quad(particle_t * particle){// -1 if no transfer, otherwise returns the rank of the reciving process
    if(particle->x > 0.5 && particle->y > 0.5){
        return (particle->quadrant_id == 3) ? -1 : 3;
    }
    else if(particle->y > 0.5){
        return (particle->quadrant_id == 2) ? -1 : 2;
    }
    else if(particle->x > 0.5){
        return (particle->quadrant_id == 1) ? -1 : 1;
    }
    return (particle->quadrant_id == 0) ? -1 : 0;
}


void free_particles(p4est_t *p4est, p4est_mesh_t * mesh) {
    for (p4est_locidx_t i = 0; i < mesh->local_num_quadrants; ++i) {
        p4est_topidx_t which_tree = -1;
        p4est_locidx_t quadrant_id = -1;
        p4est_quadrant_t *quad = p4est_mesh_quadrant_cumulative(p4est, mesh, i, &which_tree, &quadrant_id);
        particle_buffer_t *buf = (particle_buffer_t *) quad->p.user_data;
        if (buf && buf->particles) {
            free(buf->particles);
            buf->particles = NULL;
        }
    }
}

void cleanup(p4est_t *p4est) {
    p4est_destroy(p4est);
    //sc_finalize();
    sc_MPI_Finalize();
}

void print_particle_positions(p4est_t * p4est, p4est_mesh_t * mesh, mpi_context_t mpi_context, FILE *file){
    for (p4est_locidx_t i = 0; i < mesh->local_num_quadrants; ++i) {
        p4est_quadrant_t *quad = p4est_mesh_quadrant_cumulative(p4est, mesh, i, NULL, NULL);
        particle_buffer_t *buf = (particle_buffer_t *) quad->p.user_data;
        
        fprintf(file, "Rank %d: Quadrant %ld - Particles: %d\n", 
               mpi_context.mpirank, i, buf->count);
        
        // Print particle data for this quadrant
        for (int j = 0; j < buf->count; j++) {
            fprintf(file, "  Particle %d: pos=(%.6f, %.6f), vel=(%.6f, %.6f) pid: %d\n", 
                   j, buf->particles[j].x, buf->particles[j].y, 
                   buf->particles[j].vx, buf->particles[j].vy, buf->particles[j].quadrant_id);
        }
        fprintf(file, "\n");
    }
}

void add_particle(particle_t particle, particle_buffer_t *buffer){
    if(buffer == NULL || buffer->particles == NULL){
        return;
    }
    
    if(buffer->count >= buffer->capacity){
        buffer->capacity *= 2;  // Double the capacity
        buffer->particles = realloc(buffer->particles, buffer->capacity * sizeof(particle_t));
    }
    
    buffer->particles[buffer->count] = particle;
    buffer->count++;
}
void remove_particle(particle_buffer_t *buffer, int index){
    if(buffer == NULL || index < 0 || index >= buffer->count){
        return;
    }
    
    if(index < buffer->count - 1){
        buffer->particles[index] = buffer->particles[buffer->count - 1];
        memset(&buffer->particles[buffer->count - 1],/*set byte*/ 0, /*num bytes*/sizeof(particle_t));
    }
    buffer->count--;
}

//
void send_counts(particle_buffer_t ** data, MPI_Request *mpi_sends, mpi_context_t *mpi_context){
    int index = 0;
    //has access to static global_mpi_rankk
    for(int rank = 0; rank < 4; ++rank){
        if(rank == mpi_context->mpirank) continue;
        int count = (data[rank] != NULL) ? data[rank]->count : 0;
        MPI_Isend(/*Pointer to data element*/&count, 
                /*# elements*/1, 
                MPI_INT, 
                /*send to */rank, 
                /*message tag */0, 
                mpi_context->mpicomm, 
                /*promise MPI_Request*/&mpi_sends[index++]);
    }
}

void receive_counts(received_counts *counts, MPI_Request *mpi_recs, mpi_context_t *mpi_context){
    int index = 0;
    int * count_ptrs[3];
    count_ptrs[0] = &counts->count1;
    count_ptrs[1] = &counts->count2;
    count_ptrs[2] = &counts->count3;
    for(int rank = 0; rank < 4; rank++){
        if(mpi_context->mpirank == rank) continue;
        MPI_Irecv(
                count_ptrs[index],
                1,
                MPI_INT, 
                rank, 
                0,
                mpi_context->mpicomm,
                &mpi_recs[index]);
        index++;
    }

}

void loop(p4est_t *p4est, p4est_mesh_t *mesh, mpi_context_t mpi_context, int num_steps, FILE *file) {
    for (int step = 0; step < num_steps; ++step) {
        // For each local quadrant (only 1 local for -n4 demo)
        for (p4est_locidx_t i = 0; i < mesh->local_num_quadrants; ++i) {
            p4est_quadrant_t *quad = p4est_mesh_quadrant_cumulative(p4est, mesh, i, NULL, NULL);
            particle_buffer_t *buf = (particle_buffer_t *) quad->p.user_data;
            
            //4 buffers (each one for one rank)
            particle_buffer_t ** send_data = malloc(4 * sizeof(particle_buffer_t*));
            for(int j = 0; j < 4; j++){
                send_data[j] = NULL;
            }

            //initalize end pointer for particle buffer and add to conditional for loop below

            for (int j = 0; j < buf->count; ++j) {
                particle_t *p = &buf->particles[j];

                //get old distance
                double dx = p->x - planet_xyz[0];
                double dy = p->y - planet_xyz[1];
                double r = sqrt(dx*dx + dy*dy);

                //current acceleration
                double a = -grav_const * planet_mass / (r * r * r + 1e-12); //no div by 0
                double ax = a * dx;
                double ay = a * dy;

                //new velocity
                p->vx += ax * timestep;
                p->vy += ay * timestep;

                //new position
                p->x += p->vx * timestep;
                p->y += p->vy * timestep;

                //check transfer
                int new_quad = find_quad(p);
                if (new_quad != -1 && new_quad != p->quadrant_id) {
                    //initalize outgoing buffer when particle recognized
                    if(send_data[new_quad] == NULL){
                        send_data[new_quad] = malloc(sizeof(particle_buffer_t));
                        send_data[new_quad]->capacity = 10;  //Start with capacity for 10 particles
                        send_data[new_quad]->count = 0;
                        send_data[new_quad]->particles = malloc(send_data[new_quad]->capacity * sizeof(particle_t));
                    }
                    //add to outing buffer
                    add_particle(*p, send_data[new_quad]);
                    
                    //remove from self buffer
                    remove_particle(buf, j); //decremets buf->count for loop conditional
                    j--;
                }
            }
            //gone through all particles at this point, send_data buffers have been completed
            

            MPI_Request sent_promises[3];
            MPI_Request recieved_promises[3];

            received_counts rec_counts;

            receive_counts(&rec_counts, recieved_promises, &mpi_context);

            send_counts(send_data, sent_promises, &mpi_context);


            MPI_Waitall(3, recieved_promises, MPI_STATUSES_IGNORE);
            MPI_Waitall(3, sent_promises, MPI_STATUSES_IGNORE);



            //reallocate
            fprintf(file, "send_counts: 1: %d, 2: %d, 3: %d\n", rec_counts.count1, rec_counts.count2, rec_counts.count3);
            //send particles

            
            //clean up individual buffers
            for(int k = 0; k < 4; k++){
                if(send_data[k] != NULL){
                    free(send_data[k]->particles);
                    free(send_data[k]);
                }
            }
            free(send_data);
        }

        //Sync 
        MPI_Barrier(mpi_context.mpicomm);

        print_particle_positions(p4est, mesh, mpi_context, file);
    }
}

void run(int argc, char **argv) {
    mpi_context_t mpi_context = mpi_init(argc, argv);

    p4est_t *p4est = p4est_setup(&mpi_context);

    p4est_partition(p4est, 0, NULL); // 0 and NULL for uniform weight distribution

    p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    p4est_mesh_t *mesh = p4est_mesh_new(p4est, ghost, P4EST_CONNECT_FULL);

    char filename[256];
    snprintf(filename, sizeof(filename), "particles_rank%d.txt", 
             mpi_context.mpirank);

    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        exit(1);
    }
    print_particle_positions(p4est, mesh, mpi_context, file);
    //timestep = 0.1 -> 2 second demo -> 2 / 0.1 = 20
    int ns = 20;
    loop(p4est, mesh, mpi_context, ns, file);

    fclose(file);
    free_particles(p4est, mesh);
    p4est_mesh_destroy(mesh);
    p4est_ghost_destroy(ghost);

    cleanup(p4est);
}

int main(int argc, char **argv) {
    run(argc, argv);
    return 0;
}


//build data structure for sending particles out
//mpi isend number of particles to be recived to each process
//data structure for #particles notification
//mpirecive number of particles for each, keep a counter [0,3] 
//reallocate array if necessecary, (if sent particles < recived particles)
//mpi recive particle data

//work on periodic transfers later

