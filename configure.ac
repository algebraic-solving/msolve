#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([msolve], [0.1.2], [ederc@mathematik.uni-kl.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign silent-rules subdir-objects color-tests parallel-tests])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_MACRO_DIR([m4])

AM_MAINTAINER_MODE

AM_PROG_AR
LT_INIT

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([flint], [fmpz_add_ui], , [AC_MSG_ERROR([FLINT not found])])
AC_CHECK_LIB([m], [main], , [AC_MSG_ERROR([Math library not found])])
AC_CHECK_LIB([mpfr], [mpfr_init], , [AC_MSG_ERROR([MPFR not found])])
AC_CHECK_LIB([pthread], [posix_spawn])
# AC_CHECK_LIB([gomp], [omp_get_thread_num])
AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdint.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

AC_OPENMP

AX_COUNT_CPUS
AX_EXT
AX_GCC_BUILTIN

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor getdelim gettimeofday memmove memset pow sqrt strchr strstr strtol])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/fglm/Makefile
 src/neogb/Makefile
 src/usolve/Makefile
 src/msolve/Makefile
])
AC_OUTPUT
