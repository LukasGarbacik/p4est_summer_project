#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <sc_notify.h>
#include <sc_options.h>
#include <stdio.h>
#include <stdlib.h>

#include "orbit.h"

int main(int argc, char **argv) {


    printf("Hello World\n");



    return 0;
}


//main plan ->>>

// create 2d unit squre connectivity
// create mesh
// start with 4 processors for 4 squres -> then move to less processors that have more quads
// create central planet that has attractive force as particles orbit
// create random number of particels in each quadrant and fill data structures
// main loop ->>>