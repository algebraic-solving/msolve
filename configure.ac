#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([msolve], [0.9.4], [ederc@rptu.de,mohab.safey@lip6.fr])
AC_SUBST([LT_VERSION],[3:4:0])
AM_INIT_AUTOMAKE([-Wall -Werror foreign silent-rules subdir-objects color-tests parallel-tests])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_MACRO_DIR([m4])

AM_MAINTAINER_MODE

# Enable non standard exensions by defining macros such as _GNU_SOURCE
AC_USE_SYSTEM_EXTENSIONS


AM_PROG_AR
LT_INIT

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# CAUTION: Order matters for static binary building.
AC_SEARCH_LIBS([atan2],[m],
    [],
    [AC_MSG_ERROR(["The C math library was not found!"])])

AC_SEARCH_LIBS([__gmpz_init],[gmp],
    [],
    [AC_MSG_ERROR(["GMP library was not found."])])

AC_SEARCH_LIBS([mpfr_init],[mpfr],
    [],
    [AC_MSG_ERROR([MPFR library was not found.])])

AC_SEARCH_LIBS([fmpz_add_ui],[flint],
    [],
    [AC_MSG_ERROR([FLINT library was not found.])])

AC_CHECK_HEADERS([flint/nmod.h])

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
 src/fglm/Makefile
 src/neogb/Makefile
 src/usolve/Makefile
 src/msolve/Makefile
])
AC_CONFIG_LINKS([
 input_files/line_endings/in1_dos.ms:input_files/line_endings/in1_dos.ms
 input_files/line_endings/in2_dos_noeol.ms:input_files/line_endings/in2_dos_noeol.ms
 input_files/line_endings/in3_dos.ms:input_files/line_endings/in3_dos.ms
 input_files/line_endings/in4_dos.ms:input_files/line_endings/in4_dos.ms
 input_files/line_endings/in1_unix.ms:input_files/line_endings/in1_unix.ms
 input_files/line_endings/in2_unix.ms:input_files/line_endings/in2_unix.ms
 input_files/line_endings/in3_unix.ms:input_files/line_endings/in3_unix.ms
 input_files/line_endings/in4_unix.ms:input_files/line_endings/in4_unix.ms
 test/diff/diff_source.sh:test/diff/diff_source.sh
 test/diff/diff_source-gb-lm.sh:test/diff/diff_source-gb-lm.sh
 test/diff/diff_source-gb-red.sh:test/diff/diff_source-gb-red.sh
 test/diff/diff_source-gb-elim-red.sh:test/diff/diff_source-gb-elim-red.sh
 test/diff/diff_source-noarg.sh:test/diff/diff_source-noarg.sh
 test/diff/diff_source-param-fullnf-qq.sh:test/diff/diff_source-param-fullnf-qq.sh
 test/diff/diff_source-param-fullnf.sh:test/diff/diff_source-param-fullnf.sh
 test/diff/diff_source-param-nonf-qq.sh:test/diff/diff_source-param-nonf-qq.sh
 test/diff/diff_source-param-nonf.sh:test/diff/diff_source-param-nonf.sh
 test/diff/diff_source-real-param-nonf.sh:test/diff/diff_source-real-param-nonf.sh
 test/diff/diff_source-real.sh:test/diff/diff_source-real.sh
 input_files/bug-2nd-prime-bad.ms:input_files/bug-2nd-prime-bad.ms
 output_files/bug-2nd-prime-bad.g2.res:output_files/bug-2nd-prime-bad.g2.res
 test/diff/diff_bug-2nd-prime_bad.sh:test/diff/diff_bug-2nd-prime-bad.sh
 input_files/bug-68.ms:input_files/bug-68.ms
 output_files/bug-68.res:output_files/bug-68.res
 test/diff/diff_bug-68.sh:test/diff/diff_bug-68.sh
 input_files/bug-empty-tracer.ms:input_files/bug-empty-tracer.ms
 output_files/bug-empty-tracer.g2.res:output_files/bug-empty-tracer.g2.res
 test/diff/diff_bug-empty-tracer.sh:test/diff/diff_bug-empty-tracer.sh
 input_files/choice-linear-form-qq.ms:input_files/choice-linear-form-qq.ms
 output_files/choice-linear-form-qq.res:output_files/choice-linear-form-qq.res
 test/diff/diff_choice-linear-form-qq.sh:test/diff/diff_choice-linear-form-qq.sh
 input_files/choice-linear-form-qq-real-roots.ms:input_files/choice-linear-form-qq-real-roots.ms
 output_files/choice-linear-form-qq-real-roots.res:output_files/choice-linear-form-qq-real-roots.res
 test/diff/diff_choice-linear-form-qq-real-roots.sh:test/diff/diff_choice-linear-form-qq-real-roots.sh
 input_files/cp-d3-n4-p2.ms:input_files/cp-d3-n4-p2.ms
 output_files/cp-d3-n4-p2.P2.d0.res:output_files/cp-d3-n4-p2.P2.d0.res
 test/diff/diff_cp-d3-n4-p2.sh:test/diff/diff_cp-d3-n4-p2.sh
 input_files/cyclic5-16.ms:input_files/cyclic5-16.ms
 output_files/cyclic5-16.P2.d4.res:output_files/cyclic5-16.P2.d4.res
 test/diff/diff_cyclic5-16.sh:test/diff/diff_cyclic5-16.sh
 input_files/cyclic5-31.ms:input_files/cyclic5-31.ms
 output_files/cyclic5-31.P2.d4.res:output_files/cyclic5-31.P2.d4.res
 test/diff/diff_cyclic5-31.sh:test/diff/diff_cyclic5-31.sh
 input_files/cyclic5-qq.ms:input_files/cyclic5-qq.ms
 output_files/cyclic5-qq.P2.d4.res:output_files/cyclic5-qq.P2.d4.res
 output_files/cyclic5-qq.g2.res:output_files/cyclic5-qq.g2.res
 output_files/cyclic5-qq.g2.e1.res:output_files/cyclic5-qq.g2.e1.res
 output_files/cyclic5-qq.g2.e2.res:output_files/cyclic5-qq.g2.e2.res
 test/diff/diff_cyclic5-qq.sh:test/diff/diff_cyclic5-qq.sh
 input_files/eco10-31.ms:input_files/eco10-31.ms
 output_files/eco10-31.P2.d0.res:output_files/eco10-31.P2.d0.res
 test/diff/diff_eco10-31.sh:test/diff/diff_eco10-31.sh
 input_files/eco6-16.ms:input_files/eco6-16.ms
 output_files/eco6-16.P2.d4.res:output_files/eco6-16.P2.d4.res
 test/diff/diff_eco6-16.sh:test/diff/diff_eco6-16.sh
 input_files/eco6-31.ms:input_files/eco6-31.ms
 output_files/eco6-31.P2.d4.res:output_files/eco6-31.P2.d4.res
 test/diff/diff_eco6-31.sh:test/diff/diff_eco6-31.sh
 input_files/eco6-qq.ms:input_files/eco6-qq.ms
 output_files/eco6-qq.P2.d4.res:output_files/eco6-qq.P2.d4.res
 test/diff/diff_eco6-qq.sh:test/diff/diff_eco6-qq.sh
 input_files/elim-31.ms:input_files/elim-31.ms
 output_files/elim-31.g2.e1.res:output_files/elim-31.g2.e1.res
 output_files/elim-31.g2.e2.res:output_files/elim-31.g2.e2.res
 test/diff/diff_elim-31.sh:test/diff/diff_elim-31.sh
 input_files/elim-qq.ms:input_files/elim-qq.ms
 output_files/elim-qq.g2.e1.res:output_files/elim-qq.g2.e1.res
 output_files/elim-qq.g2.e2.res:output_files/elim-qq.g2.e2.res
 test/diff/diff_elim-qq.sh:test/diff/diff_elim-qq.sh
 input_files/f4sat-31.ms:input_files/f4sat-31.ms
 output_files/f4sat-31.res:output_files/f4sat-31.res
 test/diff/diff_f4sat-31.sh:test/diff/diff_f4sat-31.sh
 input_files/f4sat-byone-31.ms:input_files/f4sat-byone-31.ms
 output_files/f4sat-byone-31.res:output_files/f4sat-byone-31.res
 test/diff/diff_f4sat-byone-31.sh:test/diff/diff_f4sat-byone-31.sh
 input_files/f4sat-field-char.ms:input_files/f4sat-field-char.ms
 output_files/f4sat-field-char.res:output_files/f4sat-field-char.res
 test/diff/diff_f4sat-field-char.sh:test/diff/diff_f4sat-field-char.sh
 input_files/f4sat-is-saturated-check.ms:input_files/f4sat-is-saturated-check.ms
 output_files/f4sat-is-saturated-check.res:output_files/f4sat-is-saturated-check.res
 test/diff/diff_f4sat-is-saturated-check.sh:test/diff/diff_f4sat-is-saturated-check.sh
 input_files/f4sat-zero-input.ms:input_files/f4sat-zero-input.ms
 output_files/f4sat-zero-input.res:output_files/f4sat-zero-input.res
 test/diff/diff_f4sat-zero-input.sh:test/diff/diff_f4sat-zero-input.sh
 input_files/groebner-g2.ms:input_files/groebner-g2.ms
 output_files/groebner-g2.g2.res:output_files/groebner-g2.g2.res
 test/diff/diff_groebner-g2.sh:test/diff/diff_groebner-g2.sh
 input_files/henrion5-qq.ms:input_files/henrion5-qq.ms
 output_files/henrion5-qq.P2.d0.res:output_files/henrion5-qq.P2.d0.res
 output_files/henrion5-qq.g2.res:output_files/henrion5-qq.g2.res
 output_files/henrion5-qq.g2.e1.res:output_files/henrion5-qq.g2.e1.res
 output_files/henrion5-qq.g2.e2.res:output_files/henrion5-qq.g2.e2.res
 test/diff/diff_henrion5-qq.sh:test/diff/diff_henrion5-qq.sh
 input_files/input-overflow-16.ms:input_files/input-overflow-16.ms
 output_files/input-overflow-16.P2.d0.res:output_files/input-overflow-16.P2.d0.res
 test/diff/diff_input-overflow-16.sh:test/diff/diff_input-overflow-16.sh
 input_files/issue-230.ms:input_files/issue-230.ms
 output_files/issue-230.P1.d0.p128.res:output_files/issue-230.P1.d0.p128.res
 output_files/issue-230.P1.d0.p1024.res:output_files/issue-230.P1.d0.p1024.res
 test/diff/diff_issue-230.sh:test/diff/diff_issue-230.sh
 input_files/issue-230-squared.ms:input_files/issue-230-squared.ms
 output_files/issue-230-squared.P1.d0.p128.res:output_files/issue-230-squared.P1.d0.p128.res
 output_files/issue-230-squared.P1.d0.p1024.res:output_files/issue-230-squared.P1.d0.p1024.res
 test/diff/diff_issue-230-squared.sh:test/diff/diff_issue-230-squared.sh
 input_files/kat6-31.ms:input_files/kat6-31.ms
 output_files/kat6-31.P2.d0.res:output_files/kat6-31.P2.d0.res
 test/diff/diff_kat6-31.sh:test/diff/diff_kat6-31.sh
 input_files/kat7-qq.ms:input_files/kat7-qq.ms
 output_files/kat7-qq.P2.d0.res:output_files/kat7-qq.P2.d0.res
 test/diff/diff_kat7-qq.sh:test/diff/diff_kat7-qq.sh
 input_files/kat8-qq-truncate.ms:input_files/kat8-qq-truncate.ms
 output_files/kat8-qq-truncate.res:output_files/kat8-qq-truncate.res
 test/diff/diff_kat8-qq-truncate.sh:test/diff/diff_kat8-qq-truncate.sh
 input_files/linear-qq.ms:input_files/linear-qq.ms
 output_files/linear-qq.P2.d0.res:output_files/linear-qq.P2.d0.res
 test/diff/diff_linear-qq.sh:test/diff/diff_linear-qq.sh
 input_files/maxbitsize-bug.ms:input_files/maxbitsize-bug.ms
 output_files/maxbitsize-bug.res:output_files/maxbitsize-bug.res
 test/diff/diff_maxbitsize-bug.sh:test/diff/diff_maxbitsize-bug.sh
 input_files/mq-2-1.ms:input_files/mq-2-1.ms
 output_files/mq-2-1.res:output_files/mq-2-1.res
 test/diff/diff_mq-2-1.sh:test/diff/diff_mq-2-1.sh
 input_files/multy-16.ms:input_files/multy-16.ms
 output_files/multy-16.P2.d0.res:output_files/multy-16.P2.d0.res
 test/diff/diff_multy-16.sh:test/diff/diff_multy-16.sh
 input_files/multy-31.ms:input_files/multy-31.ms
 output_files/multy-31.P2.d0.res:output_files/multy-31.P2.d0.res
 test/diff/diff_multy-31.sh:test/diff/diff_multy-31.sh
 input_files/multy-qq.ms:input_files/multy-qq.ms
 output_files/multy-qq.P2.d0.res:output_files/multy-qq.P2.d0.res
 test/diff/diff_multy-qq.sh:test/diff/diff_multy-qq.sh
 input_files/nf-8.ms:input_files/nf-8.ms
 output_files/nf-8.res:output_files/nf-8.res
 test/diff/diff_nf-8.sh:test/diff/diff_nf-8.sh
 input_files/nf-16.ms:input_files/nf-16.ms
 output_files/nf-16.res:output_files/nf-16.res
 test/diff/diff_nf-16.sh:test/diff/diff_nf-16.sh
 input_files/nf-31.ms:input_files/nf-31.ms
 output_files/nf-31.res:output_files/nf-31.res
 test/diff/diff_nf-31.sh:test/diff/diff_nf-31.sh
 input_files/nf-lm-bug.ms:input_files/nf-lm-bug.ms
 output_files/nf-lm-bug.res:output_files/nf-lm-bug.res
 test/diff/diff_nf-lm-bug.sh:test/diff/diff_nf-lm-bug.sh
 input_files/nonradical-radicalshape-no-square-31.ms:input_files/nonradical-radicalshape-no-square-31.ms
 output_files/nonradical-radicalshape-no-square-31.c0.res:output_files/nonradical-radicalshape-no-square-31.c0.res
 test/diff/diff_nonradical-radicalshape-no-square-31.sh:test/diff/diff_nonradical-radicalshape-no-square-31.sh
 input_files/nonradical-radicalshape-no-square-qq.ms:input_files/nonradical-radicalshape-no-square-qq.ms
 output_files/nonradical-radicalshape-no-square-qq.P1.d0.p128.res:output_files/nonradical-radicalshape-no-square-qq.P1.d0.p128.res
 output_files/nonradical-radicalshape-no-square-qq.P1.d0.p1024.res:output_files/nonradical-radicalshape-no-square-qq.P1.d0.p1024.res
 test/diff/diff_nonradical-radicalshape-no-square-qq.sh:test/diff/diff_nonradical-radicalshape-no-square-qq.sh
 input_files/nonradical-radicalshape-31.ms:input_files/nonradical-radicalshape-31.ms
 output_files/nonradical-radicalshape-31.res:output_files/nonradical-radicalshape-31.res
 input_files/nonradical-radicalshape-qq.ms:input_files/nonradical-radicalshape-qq.ms
 output_files/nonradical-radicalshape-qq.P1.d0.p128.res:output_files/nonradical-radicalshape-qq.P1.d0.p128.res
 output_files/nonradical-radicalshape-qq.P1.d0.p1024.res:output_files/nonradical-radicalshape-qq.P1.d0.p1024.res
 test/diff/diff_nonradical-radicalshape-qq.sh:test/diff/diff_nonradical-radicalshape-qq.sh
 input_files/nonradical-shape-31.ms:input_files/nonradical-shape-31.ms
 output_files/nonradical-shape-31.res:output_files/nonradical-shape-31.res
 input_files/nonradical-shape-qq.ms:input_files/nonradical-shape-qq.ms
 output_files/nonradical-shape-qq.P1.d0.p128.res:output_files/nonradical-shape-qq.P1.d0.p128.res
 output_files/nonradical-shape-qq.P1.d0.p1024.res:output_files/nonradical-shape-qq.P1.d0.p1024.res
 test/diff/diff_nonradical-shape-qq.sh:test/diff/diff_nonradical-shape-qq.sh
 input_files/one-16.ms:input_files/one-16.ms
 output_files/one-16.res:output_files/one-16.res
 test/diff/diff_one-16.sh:test/diff/diff_one-16.sh
 input_files/one-31.ms:input_files/one-31.ms
 output_files/one-31.res:output_files/one-31.res
 output_files/one-31.g1.res:output_files/one-31.g1.res
 test/diff/diff_one-31.sh:test/diff/diff_one-31.sh
 input_files/one-qq.ms:input_files/one-qq.ms
 output_files/one-qq.res:output_files/one-qq.res
 output_files/one-qq.g2.res:output_files/one-qq.g2.res
 test/diff/diff_one-qq.sh:test/diff/diff_one-qq.sh
 input_files/radical-shape-31.ms:input_files/radical-shape-31.ms
 output_files/radical-shape-31.P2.d0.res:output_files/radical-shape-31.P2.d0.res
 test/diff/diff_radical-shape-31.sh:test/diff/diff_radical-shape-31.sh
 input_files/radical-shape-qq.ms:input_files/radical-shape-qq.ms
 output_files/radical-shape-qq.P2.d0.res:output_files/radical-shape-qq.P2.d0.res
 test/diff/diff_radical-shape-qq.sh:test/diff/diff_radical-shape-qq.sh
 input_files/realroot1.ms:input_files/realroot1.ms
 output_files/realroot1.p128.res:output_files/realroot1.p128.res
 output_files/realroot1.p1024.res:output_files/realroot1.p1024.res
 test/diff/diff_realroot1.sh:test/diff/diff_realroot1.sh
 input_files/realroot-extraction-exact-root.ms:input_files/realroot-extraction-exact-root.ms
 output_files/realroot-extraction-exact-root.p128.res:output_files/realroot-extraction-exact-root.p128.res
 output_files/realroot-extraction-exact-root.p1024.res:output_files/realroot-extraction-exact-root.p1024.res
 test/diff/diff_realroot-extraction-exact-root.sh:test/diff/diff_realroot-extraction-exact-root.sh
 input_files/reals-dim0-chgvar.ms:input_files/reals-dim0-chgvar.ms
 output_files/reals-dim0-chgvar.p128.res:output_files/reals-dim0-chgvar.p128.res
 output_files/reals-dim0-chgvar.p1024.res:output_files/reals-dim0-chgvar.p1024.res
 test/diff/diff_reals-dim0-chgvar.sh:test/diff/diff_reals-dim0-chgvar.sh
 input_files/reals-dim0-extract.ms:input_files/reals-dim0-extract.ms
 output_files/reals-dim0-extract.p128.res:output_files/reals-dim0-extract.p128.res
 output_files/reals-dim0-extract.p1024.res:output_files/reals-dim0-extract.p1024.res
 test/diff/diff_reals-dim0-extract.sh:test/diff/diff_reals-dim0-extract.sh
 input_files/reals-dim0.ms:input_files/reals-dim0.ms
 output_files/reals-dim0.p128.res:output_files/reals-dim0.p128.res
 output_files/reals-dim0.p1024.res:output_files/reals-dim0.p1024.res
 test/diff/diff_reals-dim0.sh:test/diff/diff_reals-dim0.sh
 input_files/reals-dim0-swapvar.ms:input_files/reals-dim0-swapvar.ms
 output_files/reals-dim0-swapvar.p128.res:output_files/reals-dim0-swapvar.p128.res
 output_files/reals-dim0-swapvar.p1024.res:output_files/reals-dim0-swapvar.p1024.res
 test/diff/diff_reals-dim0-swapvar.sh:test/diff/diff_reals-dim0-swapvar.sh
 input_files/test-lifting2.ms:input_files/test-lifting2.ms
 output_files/test-lifting2.P2.d0.res:output_files/test-lifting2.P2.d0.res
 test/diff/diff_test-lifting2.sh:test/diff/diff_test-lifting2.sh
 input_files/xy-qq.ms:input_files/xy-qq.ms
 output_files/xy-qq.P2.d0.res:output_files/xy-qq.P2.d0.res
 test/diff/diff_xy-qq.sh:test/diff/diff_xy-qq.sh
])
AC_OUTPUT
