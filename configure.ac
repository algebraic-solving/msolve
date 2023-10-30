#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([msolve], [0.6.0], [ederc@mathematik.uni-kl.de,mohab.safey@lip6.fr])
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
 msolve.pc
 include/Makefile
 src/fglm/Makefile
 src/neogb/Makefile
 src/usolve/Makefile
 src/msolve/Makefile
])
AC_CONFIG_LINKS([
 test/diff/diff_source.sh:test/diff/diff_source.sh
 input_files/cp_d_3_n_4_p_2.ms:input_files/cp_d_3_n_4_p_2.ms
 output_files/cp_d_3_n_4_p_2.res:output_files/cp_d_3_n_4_p_2.res
 test/diff/diff_cp_d_3_n_4_p_2.sh:test/diff/diff_cp_d_3_n_4_p_2.sh
 input_files/eco11-31.ms:input_files/eco11-31.ms
 output_files/eco11-31.res:output_files/eco11-31.res
 test/diff/diff_eco11-31.sh:test/diff/diff_eco11-31.sh
 input_files/elim-31.ms:input_files/elim-31.ms
 output_files/elim-31.res:output_files/elim-31.res
 output_files/elim-31.e2.res:output_files/elim-31.e2.res
 test/diff/diff_elim-31.sh:test/diff/diff_elim-31.sh
 input_files/elim-qq.ms:input_files/elim-qq.ms
 output_files/elim-qq.res:output_files/elim-qq.res
 output_files/elim-qq.e2.res:output_files/elim-qq.e2.res
 test/diff/diff_elim-qq.sh:test/diff/diff_elim-qq.sh
 input_files/F4SAT-31.ms:input_files/F4SAT-31.ms
 output_files/F4SAT-31.res:output_files/F4SAT-31.res
 test/diff/diff_F4SAT-31.sh:test/diff/diff_F4SAT-31.sh
 input_files/F4SAT-byone-31.ms:input_files/F4SAT-byone-31.ms
 output_files/F4SAT-byone-31.res:output_files/F4SAT-byone-31.res
 test/diff/diff_F4SAT-byone-31.sh:test/diff/diff_F4SAT-byone-31.sh
 input_files/kat6-31.ms:input_files/kat6-31.ms
 output_files/kat6-31.res:output_files/kat6-31.res
 test/diff/diff_kat6-31.sh:test/diff/diff_kat6-31.sh
 input_files/kat7-qq.ms:input_files/kat7-qq.ms
 output_files/kat7-qq.res:output_files/kat7-qq.res
 test/diff/diff_kat7-qq.sh:test/diff/diff_kat7-qq.sh
 input_files/multy-16.ms:input_files/multy-16.ms
 output_files/multy-16.res:output_files/multy-16.res
 test/diff/diff_multy-16.sh:test/diff/diff_multy-16.sh
 input_files/multy-31.ms:input_files/multy-31.ms
 output_files/multy-31.res:output_files/multy-31.res
 test/diff/diff_multy-31.sh:test/diff/diff_multy-31.sh
 input_files/multy-qq.ms:input_files/multy-qq.ms
 output_files/multy-qq.res:output_files/multy-qq.res
 test/diff/diff_multy-qq.sh:test/diff/diff_multy-qq.sh
 input_files/nonradical_shape-31.ms:input_files/nonradical_shape-31.ms
 output_files/nonradical_shape-31.res:output_files/nonradical_shape-31.res
 input_files/nonradical_shape-qq.ms:input_files/nonradical_shape-qq.ms
 output_files/nonradical_shape-qq.res:output_files/nonradical_shape-qq.res
 output_files/nonradical_shape-qq.p1024.res:output_files/nonradical_shape-qq.p1024.res
 test/diff/diff_nonradical_shape-qq.sh:test/diff/diff_nonradical_shape-qq.sh
 input_files/nonradical_radicalshape-31.ms:input_files/nonradical_radicalshape-31.ms
 output_files/nonradical_radicalshape-31.res:output_files/nonradical_radicalshape-31.res
 input_files/nonradical_radicalshape-qq.ms:input_files/nonradical_radicalshape-qq.ms
 output_files/nonradical_radicalshape-qq.res:output_files/nonradical_radicalshape-qq.res
 output_files/nonradical_radicalshape-qq.p3.res:output_files/nonradical_radicalshape-qq.p3.res
 test/diff/diff_nonradical_radicalshape-qq.sh:test/diff/diff_nonradical_radicalshape-qq.sh
 input_files/one-16.ms:input_files/one-16.ms
 output_files/one-16.res:output_files/one-16.res
 test/diff/diff_one-16.sh:test/diff/diff_one-16.sh
 input_files/one-31.ms:input_files/one-31.ms
 output_files/one-31.res:output_files/one-31.res
 test/diff/diff_one-31.sh:test/diff/diff_one-31.sh
 input_files/one-qq.ms:input_files/one-qq.ms
 output_files/one-qq.res:output_files/one-qq.res
 output_files/one-qq.g2.res:output_files/one-qq.g2.res
 test/diff/diff_one-qq.sh:test/diff/diff_one-qq.sh
 input_files/radical_shape-31.ms:input_files/radical_shape-31.ms
 output_files/radical_shape-31.res:output_files/radical_shape-31.res
 test/diff/diff_radical_shape-31.sh:test/diff/diff_radical_shape-31.sh
 input_files/radical_shape-qq.ms:input_files/radical_shape-qq.ms
 output_files/radical_shape-qq.res:output_files/radical_shape-qq.res
 test/diff/diff_radical_shape-qq.sh:test/diff/diff_radical_shape-qq.sh
 input_files/reals_dim0.ms:input_files/reals_dim0.ms
 output_files/reals_dim0.res:output_files/reals_dim0.res
 output_files/reals_dim0.p256.res:output_files/reals_dim0.p256.res
 test/diff/diff_reals_dim0.sh:test/diff/diff_reals_dim0.sh
 input_files/reals_dim0-swapvar.ms:input_files/reals_dim0-swapvar.ms
 output_files/reals_dim0-swapvar.res:output_files/reals_dim0-swapvar.res
 output_files/reals_dim0-swapvar.p256.res:output_files/reals_dim0-swapvar.p256.res
 test/diff/diff_reals_dim0-swapvar.sh:test/diff/diff_reals_dim0-swapvar.sh
 input_files/bug-empty-tracer.ms:input_files/bug-empty-tracer.ms
 output_files/bug-empty-tracer.res:output_files/bug-empty-tracer.res
 test/diff/diff_bug_empty_tracer.sh:test/diff/diff_bug_empty_tracer.sh
 input_files/bug-2nd-prime-bad.ms:input_files/bug-2nd-prime-bad.ms
 output_files/bug-2nd-prime-bad.res:output_files/bug-2nd-prime-bad.res
 test/diff/diff_bug_2nd_prime_bad.sh:test/diff/diff_bug_2nd_prime_bad.sh
 input_files/xy-qq.ms:input_files/xy-qq.ms
 output_files/xy-qq.res:output_files/xy-qq.res
 test/diff/diff_xy-qq.sh:test/diff/diff_xy-qq.sh
 input_files/bug_68.ms:input_files/bug_68.ms
 output_files/bug_68.res:output_files/bug_68.res
 test/diff/diff_bug_68.sh:test/diff/diff_bug_68.sh
 input_files/mq_2_1.ms:input_files/mq_2_1.ms
 output_files/mq_2_1.res:output_files/mq_2_1.res
 test/diff/diff_mq_2_1.sh:test/diff/diff_mq_2_1.sh
 input_files/input-overflow-16.ms:input_files/input-overflow-16.ms
 output_files/input-overflow-16.res:output_files/input-overflow-16.res
 test/diff/diff_input-overflow-16.sh:test/diff/diff_input-overflow-16.sh
 input_files/nf-8.ms:input_files/nf-8.ms
 output_files/nf-8.res:output_files/nf-8.res
 test/diff/diff_nf_8.sh:test/diff/diff_nf_8.sh
 input_files/nf-16.ms:input_files/nf-16.ms
 output_files/nf-16.res:output_files/nf-16.res
 test/diff/diff_nf_16.sh:test/diff/diff_nf_16.sh
 input_files/nf-31.ms:input_files/nf-31.ms
 output_files/nf-31.res:output_files/nf-31.res
 test/diff/diff_nf_31.sh:test/diff/diff_nf_31.sh
])
AC_OUTPUT
