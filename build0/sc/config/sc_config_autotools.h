/* config/sc_config_autotools.h.  Generated from sc_config_autotools.h.in by configure.  */
/* config/sc_config_autotools.h.in.  Generated from configure.ac by autoheader.  */

/* DEPRECATED (use SC_ENABLE_DEBUG instead) */
#define DEBUG 1

/* enable debug mode (assertions and extra checks) */
#define ENABLE_DEBUG 1

/* Undefine if: disable tests that use file i/o functions */
#define ENABLE_FILE_CHECKS 1

/* Undefine if: while the default alignment is sizeof (void *), this switch
   will choose the standard system malloc. For custom alignment use
   --enable-memalign=<bytes> */
#define ENABLE_MEMALIGN 1

/* Define to 1 if we are using MPI */
#define ENABLE_MPI 1

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
#define ENABLE_MPICOMMSHARED 1

/* Define to 1 if we are using MPI I/O */
#define ENABLE_MPIIO 1

/* Define to 1 if we can use MPI split nodes and shared memory */
#define ENABLE_MPISHARED 1

/* Define to 1 if we are using MPI_Init_thread */
#define ENABLE_MPITHREAD 1

/* Define to 1 if we can use MPI_Win_allocate_shared */
#define ENABLE_MPIWINSHARED 1

/* enable POSIX threads: Using --enable-pthread without arguments does not
   specify any CFLAGS; to supply CFLAGS use --enable-pthread=<PTHREAD_CFLAGS>.
   We check first for linking without any libraries and then with -lpthread;
   to avoid the latter, specify LIBS=<PTHREAD_LIBS> on configure line */
/* #undef ENABLE_PTHREAD */

/* Undefine if: disable non-thread-safe internal debug counters */
#define ENABLE_USE_COUNTERS 1

/* Undefine if: resize arrays with malloc/copy/free (HISTORIC) */
#define ENABLE_USE_REALLOC 1

/* Development with V4L2 devices works */
#define ENABLE_V4L2 1

/* Enable valgrind in executing tests */
/* #undef ENABLE_VALGRIND */

/* DEPRECATED (use SC_ENABLE_FILE_CHECKS instead) */
#define FILE_CHECKS 1

/* Define to 1 if we have MPI_Aint_diff */
#define HAVE_AINT_DIFF 1

/* Define to 1 if you have the `aligned_alloc' function. */
#define HAVE_ALIGNED_ALLOC 1

/* Define to 1 if you have the `backtrace' function. */
#define HAVE_BACKTRACE 1

/* Define to 1 if you have the `backtrace_symbols' function. */
#define HAVE_BACKTRACE_SYMBOLS 1

/* Define to 1 if you have the `basename' function. */
#define HAVE_BASENAME 1

/* Define to 1 if qsort_r conforms to BSD standard */
/* #undef HAVE_BSD_QSORT_R */

/* Define to 1 if you have the `dirname' function. */
#define HAVE_DIRNAME 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `fsync' function. */
#define HAVE_FSYNC 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if qsort_r conforms to GNU standard */
#define HAVE_GNU_QSORT_R 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if json_integer and json_real link */
/* #undef HAVE_JSON */

/* Define to 1 if you have the <libgen.h> header file. */
#define HAVE_LIBGEN_H 1

/* Define to 1 if you have the <linux/version.h> header file. */
#define HAVE_LINUX_VERSION_H 1

/* Define to 1 if you have the <linux/videodev2.h> header file. */
#define HAVE_LINUX_VIDEODEV2_H 1

/* Have we found function pthread_create. */
/* #undef HAVE_LPTHREAD */

/* Define to 1 if sqrt links successfully */
#define HAVE_MATH 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if we have MPI_INT8_T */
#define HAVE_MPI_INT8_T 1

/* Define to 1 if we have MPI_SIGNED_CHAR */
#define HAVE_MPI_SIGNED_CHAR 1

/* Define to 1 if we have MPI_UNSIGNED_LONG_LONG */
#define HAVE_MPI_UNSIGNED_LONG_LONG 1

/* Define to 1 if you have the `posix_memalign' function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the `qsort_r' function. */
#define HAVE_QSORT_R 1

/* Define to 1 if you have the <signal.h> header file. */
#define HAVE_SIGNAL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strtok_r' function. */
#define HAVE_STRTOK_R 1

/* Define to 1 if you have the `strtol' function. */
#define HAVE_STRTOL 1

/* Define to 1 if you have the `strtoll' function. */
#define HAVE_STRTOLL 1

/* Define to 1 if you have the <sys/ioctl.h> header file. */
#define HAVE_SYS_IOCTL_H 1

/* Define to 1 if you have the <sys/select.h> header file. */
#define HAVE_SYS_SELECT_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if zlib's adler32_combine links */
#define HAVE_ZLIB 1

/* minimal log priority */
/* #undef LOG_PRIORITY */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* DEPRECATED (use SC_ENABLE_MEMALIGN instead) */
#define MEMALIGN 1

/* desired alignment of allocations in bytes */
#define MEMALIGN_BYTES (8)

/* DEPRECATED (use SC_ENABLE_MPI instead) */
#define MPI 1

/* DEPRECATED (use SC_ENABLE_MPIIO instead) */
#define MPIIO 1

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Name of package */
#define PACKAGE "libsc"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "p4est@ins.uni-bonn.de"

/* Define to the full name of this package. */
#define PACKAGE_NAME "libsc"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "libsc 2.8.7"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "libsc"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.8.7"

/* DEPRECATED (use SC_WITH_PAPI instead) */
/* #undef PAPI */

/* Use builtin getopt */
/* #undef PROVIDE_GETOPT */

/* DEPRECATED (use SC_ENABLE_PTHREAD instead) */
/* #undef PTHREAD */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* DEPRECATED (use SC_ENABLE_USE_COUNTERS instead) */
#define USE_COUNTERS 1

/* DEPRECATED (use SC_ENABLE_USE_REALLOC instead) */
#define USE_REALLOC 1

/* Define to 1 if using autoconf build */
#define USING_AUTOCONF 1

/* Version number of package */
#define VERSION "2.8.7"

/* Package major version */
#define VERSION_MAJOR 2

/* Package minor version */
#define VERSION_MINOR 8

/* Package point version */
#define VERSION_POINT 7

/* enable Flop counting with papi */
/* #undef WITH_PAPI */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */
