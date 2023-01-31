#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([msolve], [0.4.9], [ederc@mathematik.uni-kl.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign silent-rules subdir-objects color-tests parallel-tests])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_MACRO_DIR([m4])

AM_MAINTAINER_MODE

AM_PROG_AR
LT_INIT

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# CAUTION: Order matters for static binary building.
AC_CHECK_LIB([m], [main], , [AC_MSG_ERROR([Math library not found])])
AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])
AC_CHECK_LIB([mpfr], [mpfr_init], , [AC_MSG_ERROR([MPFR not found])])
AC_CHECK_LIB([flint], [fmpz_add_ui], , [AC_MSG_ERROR([FLINT not found])])
# AC_CHECK_LIB([pthread], [posix_spawn])

# check if we want OpenMP support
AC_ARG_ENABLE([openmp],
	[  --enable-openmp         Enable OpenMP support],
	[case "${enableval}" in
		yes) 	openmp=true ;;
		no)		openmp=false ;;
		*)		AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
	esac],[openmp=true])

if test x$openmp = xtrue ; then
		AC_OPENMP
fi

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


AX_COUNT_CPUS
AX_EXT
AX_GCC_BUILTIN([__builtin_constant_p])
AX_GCC_BUILTIN([__builtin_clzll])
AX_GCC_BUILTIN([__builtin_clzl])

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
