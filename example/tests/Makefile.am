# Makefile.am in example/tests
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_3D

## Define the program to build. It will be installed in the bin directory.
bin_PROGRAMS += example/tests/working

## Define the source file for the 'working' program.
example_tests_working_SOURCES = example/tests/working.c

endif

## does not require the 3d stuff at all
bin_PROGRAMS += example/tests/orbit
example_tests_orbit_SOURCES = example/tests/orbit.c