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

#define p_per_quad 5

static sc_rand_state_t global_rand_state;

mpi_context_t mpi_init(int argc, char **argv) {

    mpi_context_t mpi_context;

    int mpiret = sc_MPI_Init(&argc, &argv);

    mpi_context.mpicomm = sc_MPI_COMM_WORLD;
    sc_MPI_Comm_rank(mpi_context.mpicomm, &mpi_context.mpirank);
    sc_MPI_Comm_size(mpi_context.mpicomm, &mpi_context.mpisize);

    SC_CHECK_MPI (mpiret);

    sc_init (mpi_context.mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

    global_rand_state = (uint64_t)(time(NULL) + mpi_context.mpirank);

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
    
    // Add some randomness within the quadrant
    double rand_x = sc_rand(&global_rand_state) / 2;
    double rand_y = sc_rand(&global_rand_state) / 2;
    
    particle.x = xs + rand_x;
    particle.y = ys + rand_y;
    
    //velocity of particle (0 temp)
    particle.vx = 0;
    particle.vy = 0;
    
    return particle;
}

void quad_init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant){
    particle_buffer_t *buf = (particle_buffer_t *) quadrant->p.user_data;
    buf->count = p_per_quad;
    buf->particles = malloc(buf->count * sizeof(particle_t));
    for(int i = 0; i < p_per_quad; i++){
        buf->particles[i] = particle_single_init(p4est, which_tree, quadrant);
    }
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
    sc_finalize();
    sc_MPI_Finalize();
}

void print_particle_positions(p4est_t * p4est, p4est_mesh_t * mesh, mpi_context_t mpi_context){
    for (p4est_locidx_t i = 0; i < mesh->local_num_quadrants; ++i) {
        p4est_topidx_t which_tree = -1;
        p4est_locidx_t quadrant_id = -1;
        p4est_quadrant_t *quad = p4est_mesh_quadrant_cumulative(p4est, mesh, i, &which_tree, &quadrant_id);
        particle_buffer_t *buf = (particle_buffer_t *) quad->p.user_data;

        char filename[256];
        snprintf(filename, sizeof(filename), "particles_rank%d_quadrant%d.txt", mpi_context.mpirank, i);

        FILE *file = fopen(filename, "w");
        if (!file) {
            printf("Error: Could not open file %s\n", filename);
            continue;
        }
        
        fprintf(file, "Rank %d: Quadrant %ld - Particles: %d\n", 
               mpi_context.mpirank, i, buf->count);
        
        // Print particle data for this quadrant
        for (int j = 0; j < buf->count; j++) {
            fprintf(file, "  Particle %d: pos=(%.6f, %.6f), vel=(%.6f, %.6f)\n", 
                   j, buf->particles[j].x, buf->particles[j].y, 
                   buf->particles[j].vx, buf->particles[j].vy);
        }
        fclose(file);
    }
}



void run(int argc, char **argv) {
    mpi_context_t mpi_context = mpi_init(argc, argv);

    p4est_t *p4est = p4est_setup(&mpi_context);

    p4est_partition(p4est, 0, NULL); // 0 and NULL for uniform weight distribution

    p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    p4est_mesh_t *mesh = p4est_mesh_new(p4est, ghost, P4EST_CONNECT_FULL);

    print_particle_positions(p4est, mesh, mpi_context);

    

    free_particles(p4est, mesh);
    p4est_mesh_destroy(mesh);
    p4est_ghost_destroy(ghost);

    cleanup(p4est);

}
int main(int argc, char **argv) {
    run(argc, argv);
    return 0;
}


//next
//define central planet
//define particle velocities
//main timestep loop
