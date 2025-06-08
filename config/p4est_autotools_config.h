/* config/p4est_autotools_config.h.  Generated from p4est_autotools_config.h.in by configure.  */
/* config/p4est_autotools_config.h.in.  Generated from configure.ac by autoheader.  */

/* DEPRECATED (use P4EST_ENABLE_BUILD_2D instead) */
#define BUILD_2D 1

/* DEPRECATED (use P4EST_ENABLE_BUILD_3D instead) */
#define BUILD_3D 1

/* DEPRECATED (use P4EST_ENABLE_BUILD_P6EST instead) */
#define BUILD_P6EST 1

/* DEPRECATED (use P4EST_ENABLE_DEBUG instead) */
/* #undef DEBUG */

/* Undefine if: disable the 2D library */
#define ENABLE_BUILD_2D 1

/* Undefine if: disable the 3D library */
#define ENABLE_BUILD_3D 1

/* Undefine if: disable hybrid 2D+1D p6est library */
#define ENABLE_BUILD_P6EST 1

/* enable debug mode (assertions and extra checks) */
/* #undef ENABLE_DEBUG */

/* Undefine if: disable tests that use file i/o functions */
#define ENABLE_FILE_CHECKS 1

/* use depreacted data file format */
/* #undef ENABLE_FILE_DEPRECATED */

/* Undefine if: while the default alignment is sizeof (void *), this switch
   will choose the standard system malloc. For custom alignment use
   --enable-memalign=<bytes> */
#define ENABLE_MEMALIGN 1

/* Define to 1 if we are using MPI */
/* #undef ENABLE_MPI */

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
/* #undef ENABLE_MPICOMMSHARED */

/* Define to 1 if we are using MPI I/O */
/* #undef ENABLE_MPIIO */

/* Define to 1 if we can use MPI split nodes and shared memory */
/* #undef ENABLE_MPISHARED */

/* Define to 1 if we are using MPI_Init_thread */
/* #undef ENABLE_MPITHREAD */

/* Define to 1 if we can use MPI_Win_allocate_shared */
/* #undef ENABLE_MPIWINSHARED */

/* enable POSIX threads: Using --enable-pthread without arguments does not
   specify any CFLAGS; to supply CFLAGS use --enable-pthread=<PTHREAD_CFLAGS>.
   We check first for linking without any libraries and then with -lpthread;
   to avoid the latter, specify LIBS=<PTHREAD_LIBS> on configure line */
/* #undef ENABLE_PTHREAD */

/* Development with V4L2 devices works */
/* #undef ENABLE_V4L2 */

/* Enable valgrind in executing tests */
/* #undef ENABLE_VALGRIND */

/* Undefine if: write vtk ascii file data */
#define ENABLE_VTK_BINARY 1

/* Undefine if: disable zlib compression for vtk binary data */
#define ENABLE_VTK_COMPRESSION 1

/* use doubles for vtk file data */
/* #undef ENABLE_VTK_DOUBLES */

/* DEPRECATED (use P4EST_ENABLE_FILE_CHECKS instead) */
#define FILE_CHECKS 1

/* DEPRECATED (use P4EST_ENABLE_FILE_DEPRECATED instead) */
/* #undef FILE_DEPRECATED */

/* Define to 1 if we have MPI_Aint_diff */
/* #undef HAVE_AINT_DIFF */

/* Define to 1 if you have the `aligned_alloc' function. */
#define HAVE_ALIGNED_ALLOC 1

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* Define to 1 if qsort_r conforms to BSD standard */
/* #undef HAVE_BSD_QSORT_R */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if qsort_r conforms to GNU standard */
#define HAVE_GNU_QSORT_R 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if json_integer and json_real link */
/* #undef HAVE_JSON */

/* Have we found function pthread_create. */
/* #undef HAVE_LPTHREAD */

/* Define to 1 if sqrt links successfully */
#define HAVE_MATH 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if we have MPI_INT8_T */
/* #undef HAVE_MPI_INT8_T */

/* Define to 1 if we have MPI_SIGNED_CHAR */
/* #undef HAVE_MPI_SIGNED_CHAR */

/* Define to 1 if we have MPI_UNSIGNED_LONG_LONG */
/* #undef HAVE_MPI_UNSIGNED_LONG_LONG */

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1

/* Define to 1 if you have the `posix_memalign' function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if zlib's adler32_combine links */
#define HAVE_ZLIB 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* DEPRECATED (use P4EST_ENABLE_MEMALIGN instead) */
#define MEMALIGN 1

/* desired alignment of allocations in bytes */
#define MEMALIGN_BYTES (8)

/* DEPRECATED (use P4EST_WITH_METIS instead) */
/* #undef METIS */

/* DEPRECATED (use P4EST_ENABLE_MPI instead) */
/* #undef MPI */

/* DEPRECATED (use P4EST_ENABLE_MPIIO instead) */
/* #undef MPIIO */

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Name of package */
#define PACKAGE "p4est"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "p4est@ins.uni-bonn.de"

/* Define to the full name of this package. */
#define PACKAGE_NAME "p4est"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "p4est 2.8.7.12-cf1ab"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "p4est"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.8.7.12-cf1ab"

/* DEPRECATED (use P4EST_WITH_PETSC instead) */
/* #undef PETSC */

/* Use builtin getopt */
/* #undef PROVIDE_GETOPT */

/* DEPRECATED (use P4EST_ENABLE_PTHREAD instead) */
/* #undef PTHREAD */

/* DEPRECATED (use P4EST_WITH_SC instead) */
/* #undef SC */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if using autoconf build */
#define USING_AUTOCONF 1

/* Version number of package */
#define VERSION "2.8.7.12-cf1ab"

/* Package major version */
#define VERSION_MAJOR 2

/* Package minor version */
#define VERSION_MINOR 8

/* Package point version */
#define VERSION_POINT 7.12-cf1ab

/* DEPRECATED (use P4EST_ENABLE_VTK_BINARY instead) */
#define VTK_BINARY 1

/* DEPRECATED (use P4EST_ENABLE_VTK_COMPRESSION instead) */
#define VTK_COMPRESSION 1

/* DEPRECATED (use P4EST_ENABLE_VTK_DOUBLES instead) */
/* #undef VTK_DOUBLES */

/* enable metis-dependent code */
/* #undef WITH_METIS */

/* enable PETSc-dependent code */
/* #undef WITH_PETSC */

/* specifiy how to depend on package sc (optional). If the option value is
   literal no or the option is not present, use the source subdirectory. If
   the option value is the literal external, expect all necessary environment
   variables to be set to compile and link against sc and to run the examples.
   Otherwise, the option value must be an absolute path to the toplevel
   directory of a make installed sc. */
/* #undef WITH_SC */
