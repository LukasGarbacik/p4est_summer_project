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
#include <stdbool.h>
#include <math.h>

#include "orbit.h"

#define p_per_quad 5

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

/*void check_and_fix_periodic(particle_t * particle){
    bool fixed = false;
    if (particle->x < 0 || particle->x >= 1 || particle->y < 0 || particle->y >= 1) {
        fixed = true;
        //180 degree rotation formula
        particle->x = 1.0 - particle->x;
        particle->y = 1.0 - particle->y;
        //fix domin to unit square
        if (particle->x < 0) particle->x += 1.0;
        if (particle->x >= 1) particle->x -= 1.0;
        if (particle->y < 0) particle->y += 1.0;
        if (particle->y >= 1) particle->y -= 1.0;
    }
    if (fixed) {
        particle->vx = -particle->vx;
        particle->vy = -particle->vy;
    }
}*/
void check_and_fix_periodic(particle_t * particle){
    if(particle->x < 0){
        particle->x = 1 - fabs(particle->x);
        particle->y = 1 - particle->y;

        particle->vy = -particle->vy; //vy flipped
    }
    else if(particle->x > 1){
        particle->x = particle->x - 1;
        particle->y = 1 - particle->y;

        particle->vy = -particle->vy; //vy flipped
    }
    else if(particle->y < 0){
        particle->x = 1 - particle->x;
        particle->y = 1 - fabs(particle->y);

        particle->vx = -particle->vx; //vx flipped
    }
    else if(particle->y > 1){
        particle->x = 1 - particle->x;
        particle->y = particle->y - 1;

        particle->vx = -particle->vx; //vx flipped
    }
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
    //consider periodic movement and correct if necessary
    check_and_fix_periodic(particle);

    int quad = 0;
    if (particle->x > 0.5 && particle->y > 0.5) quad = 3;
    else if (particle->y > 0.5) quad = 2;
    else if (particle->x > 0.5) quad = 1;
    else quad = 0;

    return (quad == particle->quadrant_id) ? -1 : quad;
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

void print_particle_positions(p4est_t * p4est, p4est_mesh_t * mesh, mpi_context_t mpi_context, FILE *file, int ln){
    for (p4est_locidx_t i = 0; i < mesh->local_num_quadrants; ++i) {
        p4est_quadrant_t *quad = p4est_mesh_quadrant_cumulative(p4est, mesh, i, NULL, NULL);
        particle_buffer_t *buf = (particle_buffer_t *) quad->p.user_data;
        
        fprintf(file, "Rank %d: Quadrant %ld - Particles: %d\n", 
               mpi_context.mpirank, i, buf->count);
        
        // Print particle data for this quadrant
        for (int j = 0; j < buf->count; j++) {
            fprintf(file, "  Particle %d: pos=(%.6f, %.6f), vel=(%.6f, %.6f) loop #: %d\n", 
                   j, buf->particles[j].x, buf->particles[j].y, 
                   buf->particles[j].vx, buf->particles[j].vy, ln);
        }
        fprintf(file, "\n");
    }
}

void add_particle(particle_t particle, particle_buffer_t *buffer){
    if(buffer == NULL || buffer->particles == NULL){
        return;
    }
    if(global_mpi_rank == 0){
        printf("this is where the particle is getting added to an outgoing buffer. self-rank: %d\n\n", global_mpi_rank);
    }
    if(buffer->count >= buffer->capacity){
        buffer->capacity *= 2;  // Double the capacity
        buffer->particles = realloc(buffer->particles, buffer->capacity * sizeof(particle_t));
    }
    
    buffer->particles[buffer->count] = particle;
    buffer->count++;
}
void remove_particle(particle_buffer_t *buffer, int index, int ln){
    if(buffer == NULL || index < 0 || index >= buffer->count){
        return;
    }
    if(global_mpi_rank == 0){
        printf("particle is getting removed from outgoing ranks self buffer. self-rank: %d\n\n", global_mpi_rank);
        printf("particle pos(%3f, %3f) pid: %d, loop number: %d\n\n", buffer->particles[index].x, buffer->particles[index].y, buffer->particles[index].quadrant_id, ln);
        printf("buffer_count: %d, arg index: %d\n", buffer->count, index);
        printf("HERE^^\n");
    }
    
    if(index < buffer->count){
        if(global_mpi_rank == 0) printf("GOT HERE\n");
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
        if(count != 0 && global_mpi_rank == 0){
            printf("COUNT OF PARTICLES SENT TO RANK: %d = %d", rank, count);
        }
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
    count_ptrs[0] = &counts->count1[0];
    count_ptrs[1] = &counts->count2[0];
    count_ptrs[2] = &counts->count3[0];
    for(int rank = 0; rank < 4; rank++){
        if(mpi_context->mpirank == rank) continue;
        count_ptrs[index][1] = rank; //save rank to order particle placement in buffer
        MPI_Irecv(
                count_ptrs[index],  //pass the pointer
                1,
                MPI_INT, 
                /* receive rank*/rank, 
                0,
                mpi_context->mpicomm,
                &mpi_recs[index]);
        index++;
    }

}

void reallocate_buffer(particle_buffer_t * buffer, received_counts * counts){
    int new_total = buffer->count + counts->count1[0] + counts->count2[0] + counts->count3[0];
    if(global_mpi_rank == 0) printf("reallocation: new_total: %d, buf cap: %d", new_total, buffer->capacity);
    if(new_total > buffer->capacity){
        buffer->capacity *= 2;
        buffer->particles = realloc(buffer->particles, buffer->capacity * sizeof(particle_t));
    }
    else if(new_total < buffer->capacity / 4 && new_total > 10){ //only deallocate for higher number of particles
        buffer->capacity /= 2;
        buffer->particles = realloc(buffer->particles, buffer->capacity * sizeof(particle_t));
    }
}

int get_prefix_helper(particle_buffer_t * buffer, received_counts * counts, int rank){//index place into buffer
    int ret = (buffer != NULL) ? buffer->count : 0;
    //reminder: count[1-3] holds the particle count inc. from corresponding rank: <count, inc. rank> 
    if(counts->count1[1] == rank){
        return ret;
    }
    else if(counts->count2[1] == rank){
        return ret + counts->count1[0];
    }
    return ret + counts->count1[0] + counts->count2[0];
}

void receive_particles(particle_buffer_t * new_buffer,  MPI_Request *mpi_recs, received_counts *counts ,  mpi_context_t * mpi_context){
    //loop through ranks, match
    int index = 0;
    for(int rank = 0; rank < 4; rank++){
        if(rank == mpi_context->mpirank) continue;

        // Find the count for this rank
        int count = 0;
        if (counts->count1[1] == rank) count = counts->count1[0];
        else if (counts->count2[1] == rank) count = counts->count2[0];
        else if (counts->count3[1] == rank) count = counts->count3[0];

        if(global_mpi_rank == 1 && count > 0) printf("rank 1 receiving particle here");

        if (count > 0) {
            int prefix = get_prefix_helper(new_buffer, counts, rank);
            MPI_Irecv(
                &new_buffer->particles[prefix], //destination
                count * sizeof(particle_t), //number of bytes
                MPI_BYTE,
                rank,
                1, //tag
                mpi_context->mpicomm,
                &mpi_recs[index]
            );
        } else {
            // If no particles expected, set request to MPI_REQUEST_NULL
            mpi_recs[index] = MPI_REQUEST_NULL;
        }
        index++;
    }
}

void send_particles(particle_buffer_t ** send_data, MPI_Request *mpi_sends, mpi_context_t *mpi_context){
    int index = 0;
    for(int rank = 0; rank < 4; rank++){
        if(rank == mpi_context->mpirank) continue;
        
        if(send_data[rank] != NULL && send_data[rank]->count > 0) {
            MPI_Isend(
                send_data[rank]->particles,
                send_data[rank]->count * sizeof(particle_t),  //byte size of outgoing particle buffer
                MPI_BYTE, 
                rank, 
                1,  //same tag as receieve
                mpi_context->mpicomm, 
                &mpi_sends[index]
            );
        } else {
            mpi_sends[index] = MPI_REQUEST_NULL;//null request to pass checker
        }
        index++;
    }
}



void loop(p4est_t *p4est, p4est_mesh_t *mesh, mpi_context_t mpi_context, int num_steps, FILE *file) {
    for (int step = 0; step < num_steps; ++step) {
        // For each local quadrant (only 1 local for -n4 demo)
        for (p4est_locidx_t i = 0; i < mesh->local_num_quadrants; ++i) {
            p4est_quadrant_t *quad = p4est_mesh_quadrant_cumulative(p4est, mesh, i, NULL, NULL);
            particle_buffer_t *buf = (particle_buffer_t *) quad->p.user_data;

            bool particle_was_transfered = false;
            
            //4 buffers (each one for one rank)
            particle_buffer_t ** send_data = malloc(4 * sizeof(particle_buffer_t*));
            for(int j = 0; j < 4; j++){
                send_data[j] = NULL;
            }

            //initalize end pointer for particle buffer and add to conditional for loop below
            int DEBUG_PARTICLE_TRANSFER_COUNT = 0;
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
                    DEBUG_PARTICLE_TRANSFER_COUNT++;

                    if (mpi_context.mpirank == 0) {
                        printf("[SEND] Rank 0 sending particle to rank %d: pos=(%.6f, %.6f), vel=(%.6f, %.6f), pid=%d, loop=%d\n",
                            new_quad, p->x, p->y, p->vx, p->vy, p->quadrant_id, step);
                    }
                    p->quadrant_id = new_quad;
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
                    remove_particle(buf, j, step); //decremets buf->count for loop conditional
                    j--;
                }
            }
            //gone through all particles at this point, send_data buffers have been completed
            
            //do count transfer no matter what
            MPI_Request sent_promises[3];
            MPI_Request recieved_promises[3];

            received_counts rec_counts;

            receive_counts(&rec_counts, recieved_promises, &mpi_context);

            send_counts(send_data, sent_promises, &mpi_context);

            fprintf(file, "particles exiting rank: %d = %d", mpi_context.mpirank, DEBUG_PARTICLE_TRANSFER_COUNT);

            MPI_Waitall(3, recieved_promises, MPI_STATUSES_IGNORE);
            MPI_Waitall(3, sent_promises, MPI_STATUSES_IGNORE);
            if(rec_counts.count1[0] + rec_counts.count2[0] + rec_counts.count3[0] > 0 && global_mpi_rank == 1){
                printf("Rank %d received counts: [%d from %d] [%d from %d] [%d from %d]\n",
                    mpi_context.mpirank,
                    rec_counts.count1[0], rec_counts.count1[1],
                    rec_counts.count2[0], rec_counts.count2[1],
                    rec_counts.count3[0], rec_counts.count3[1]);
            }
            
            particle_was_transfered = (rec_counts.count1[0] + rec_counts.count2[0] + rec_counts.count3[0] > 0) ? true : false;

            //always do particle communicatioin to stay in sync
            MPI_Request particle_r_promises[3];
            MPI_Request particle_s_promises[3];
            
            if(particle_was_transfered){
                reallocate_buffer(buf, &rec_counts);
            }
            
            receive_particles(buf, particle_r_promises, &rec_counts, &mpi_context);
            send_particles(send_data, particle_s_promises, &mpi_context);

            MPI_Waitall(3, particle_r_promises, MPI_STATUSES_IGNORE);
            MPI_Waitall(3, particle_s_promises, MPI_STATUSES_IGNORE);
            
            if(particle_was_transfered){
                buf->count += rec_counts.count1[0] + rec_counts.count2[0] + rec_counts.count3[0];
            }

            // Print received particles for rank 1
            if (mpi_context.mpirank == 1 && particle_was_transfered) {
                for (int k = 0; k < buf->count; ++k) {
                    particle_t *recv_p = &buf->particles[k];
                    printf("[RECV] Rank 1 new buffer particle (%d): pos=(%.6f, %.6f), vel=(%.6f, %.6f), pid=%d, loop=%d\n",
                        k, recv_p->x, recv_p->y, recv_p->vx, recv_p->vy, recv_p->quadrant_id, step);
                    printf("buffer info after send// count: %d, capacity: %d", buf->count, buf->capacity);
                }
            }
            
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

        print_particle_positions(p4est, mesh, mpi_context, file, step);
    }
}

void run(int argc, char **argv) {
    mpi_context_t mpi_context = mpi_init(argc, argv);

    p4est_t *p4est = p4est_setup(&mpi_context);

    p4est_partition(p4est, 0, NULL); // 0 and NULL for uniform weight distribution

    p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    p4est_mesh_t *mesh = p4est_mesh_new(p4est, ghost, P4EST_CONNECT_FULL);

    //GDB hold start for debugging (new terminal) gdb ./orbit PID
    
    /*if (mpi_context.mpirank == 0) { // or just always, if not using MPI
        printf("PID %d: Press Enter to continue...\n", getpid());
        fflush(stdout);
        getchar();
    }*/

    char filename[256];
    snprintf(filename, sizeof(filename), "particles_rank%d.txt", 
             mpi_context.mpirank);

    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        exit(1);
    }
    print_particle_positions(p4est, mesh, mpi_context, file, 0);
    //timestep = 0.1 -> 2 second demo -> 2 / 0.1 = 20
    int ns = 200;
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

