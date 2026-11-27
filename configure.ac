#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([msolve], [0.11.0], [ederc@rptu.de,mohab.safey@lip6.fr])
AC_SUBST([LT_VERSION],[4:0:0])
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
AC_CHECK_HEADERS([inttypes.h stdint.h])

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
# Record the SIMD kernels enabled for this compilation in the configure summary.
msolve_simd_kernels=
AC_DEFUN([MSOLVE_CHECK_TARGET_FEATURE], [
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#ifndef $1
#error Target feature is not enabled
#endif
]], [[]])], [
    m4_ifval([$3], [AC_DEFINE([$3], [1], [Enable $2 code])])
    msolve_simd_kernels="$msolve_simd_kernels $2"
  ])
])

# MSOLVE_AVX_OPTION(option, name, config macro, target condition, intrinsic test,
#                  requirements for an explicit enable)
AC_DEFUN([MSOLVE_AVX_OPTION], [
  AC_ARG_ENABLE([$1],
      [AS_HELP_STRING([--enable-$1], [enable $2 kernels (default: compiler target)])],
      [], [enable_$1=auto])
  AS_CASE([$enable_$1], [yes|no|auto], [],
      [AC_MSG_ERROR([bad value $enable_$1 for --enable-$1])])
  AS_IF([test "x$enable_$1" != xno], [
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if !($4)
#error Required kernel features are not enabled
#endif
#include <immintrin.h>
]], [[$5]])], [
      AC_DEFINE([$3], [1], [Enable $2 code])
      msolve_simd_kernels="$msolve_simd_kernels $2"
      enable_$1=yes
    ], [
      AS_IF([test "x$enable_$1" = xyes], [AC_MSG_ERROR([--enable-$1 requires $6])])
      enable_$1=no
    ])
  ])
])

MSOLVE_CHECK_TARGET_FEATURE([__aarch64__], [NEON])
MSOLVE_CHECK_TARGET_FEATURE([__SSE2__], [SSE2], [HAVE_SSE2])
MSOLVE_CHECK_TARGET_FEATURE([__SSSE3__], [SSSE3], [HAVE_SSSE3])
MSOLVE_AVX_OPTION([avx2], [AVX2], [HAVE_AVX2], [defined(__AVX2__)], [
    __m256i v = _mm256_set1_epi32(1);
    return _mm256_extract_epi32(_mm256_add_epi32(v, v), 0);
], [an AVX2 compiler target (for example, -mavx2 in CFLAGS)])
AS_IF([test "x$enable_avx2" = xyes], [
  MSOLVE_CHECK_TARGET_FEATURE([__FMA__], [FMA3], [HAVE_FMA3])
])
MSOLVE_AVX_OPTION([avx512f], [AVX-512F], [HAVE_AVX512_F], [defined(__AVX512F__)], [
    __m512i v = _mm512_set1_epi32(1);
    return _mm512_cvtsi512_si32(_mm512_add_epi32(v, v));
], [an AVX-512F compiler target (for example, -mavx512f in CFLAGS)])
MSOLVE_AVX_OPTION([avx512bw], [AVX-512BW], [HAVE_AVX512_BW],
    [defined(HAVE_AVX512_F) && defined(__AVX512BW__)], [
    __m512i v = _mm512_set1_epi16(1);
    return _mm512_cvtsi512_si32(_mm512_mullo_epi16(v, v));
], [enabled AVX-512F kernels and an AVX-512BW compiler target (for example, -mavx512bw in CFLAGS)])
AS_IF([test "x$enable_avx512f" = xyes], [
  MSOLVE_CHECK_TARGET_FEATURE([__AVX512DQ__], [AVX-512DQ], [HAVE_AVX512_DQ])
  MSOLVE_CHECK_TARGET_FEATURE([__AVX512IFMA__], [AVX-512IFMA], [HAVE_AVX512_IFMA])
], [test "x$enable_avx2" = xyes], [
  MSOLVE_CHECK_TARGET_FEATURE([__AVXIFMA__], [AVX-IFMA])
])
AX_GCC_BUILTIN([__builtin_constant_p])
AX_GCC_BUILTIN([__builtin_clzll])
AX_GCC_BUILTIN([__builtin_clzl])

# FLINT's build configuration does not guarantee that its fft_small headers
# are usable with our compiler target (for example, AVX2 also needs FMA).
# Compile the vector operations used by Hensel lifting without adding flags.
AC_CACHE_CHECK([whether FLINT fft_small is usable with the compiler target],
    [msolve_cv_flint_fft_small], [
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <flint/flint.h>
#if __FLINT_RELEASE < 30200 || !defined(FLINT_HAVE_FFT_SMALL) \
    || !FLINT_HAVE_FFT_SMALL
#error FLINT fft_small is unavailable
#else
#include <flint/fft_small.h>

/* External inputs and an observable result prevent constant folding from
   hiding unsupported intrinsics in the inline vector operations. */
void msolve_check_fft_small(double *out, const double *a, const double *b,
                           double p, double pinv)
{
    vec4d va = vec4d_load_aligned(a), vb = vec4d_load_aligned(b);
    vec4d vp = vec4d_set_d(p), vi = vec4d_set_d(pinv);
    vec4d h = vec4d_mul(va, vb), l = vec4d_fmsub(va, vb, h);
    vec4d r = vec4d_add(vec4d_fnmadd(vec4d_round(vec4d_mul(h, vi)), vp, h), l);
    vec4d_store_aligned(out,
        vec4d_fnmadd(vec4d_round(vec4d_mul(r, vi)), vp, r));
}
#endif
  ]])], [msolve_cv_flint_fft_small=yes], [msolve_cv_flint_fft_small=no])
])
AS_IF([test "x$msolve_cv_flint_fft_small" = xyes], [
  AC_DEFINE([HAVE_FLINT_FFT_SMALL], [1],
      [Define if FLINT fft_small is usable with the compiler target])
])

# Checks for library functions.
AC_CHECK_FUNCS([floor getdelim gettimeofday memmove memset malloc realloc pow sqrt strchr strstr strtol])

# Configure version numbers
m4_define([version_split], m4_split(AC_PACKAGE_VERSION, [[.]]))
m4_define([version_major], m4_argn(1, version_split))
m4_define([version_minor], m4_argn(2, version_split))
m4_define([version_patch], m4_argn(3, version_split))

AC_SUBST([MSOLVE_VERSION_MAJOR],[version_major])
AC_SUBST([MSOLVE_VERSION_MINOR],[version_minor])
AC_SUBST([MSOLVE_VERSION_PATCH],[version_patch])


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 msolve.pc
 src/fglm/Makefile
 src/neogb/Makefile
 src/usolve/Makefile
 src/msolve/Makefile
 src/msolve/msolve-version.h
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
 test/diff/diff_line-endings.sh:test/diff/diff_line-endings.sh
 test/diff/diff_source.sh:test/diff/diff_source.sh
 test/diff/diff_source-gb-lm.sh:test/diff/diff_source-gb-lm.sh
 test/diff/diff_source-gb-red.sh:test/diff/diff_source-gb-red.sh
 test/diff/diff_source-gb-elim-drl.sh:test/diff/diff_source-gb-elim-drl.sh
 test/diff/diff_source-gb-elim-red.sh:test/diff/diff_source-gb-elim-red.sh
 test/diff/diff_source-gb-lex-lm.sh:test/diff/diff_source-gb-lex-lm.sh
 test/diff/diff_source-gb-lex-red.sh:test/diff/diff_source-gb-lex-red.sh
 test/diff/diff_source-gb-lex-red-qq.sh:test/diff/diff_source-gb-lex-red-qq.sh
 test/diff/diff_source-noarg.sh:test/diff/diff_source-noarg.sh
 test/diff/diff_source-param-fullnf-qq.sh:test/diff/diff_source-param-fullnf-qq.sh
 test/diff/diff_source-param-fullnf.sh:test/diff/diff_source-param-fullnf.sh
 test/diff/diff_source-param-nonf-qq.sh:test/diff/diff_source-param-nonf-qq.sh
 test/diff/diff_source-param-nonf.sh:test/diff/diff_source-param-nonf.sh
 test/diff/diff_source-real-param-nonf.sh:test/diff/diff_source-real-param-nonf.sh
 test/diff/diff_source-real.sh:test/diff/diff_source-real.sh
 input_files/boolean-elim10.ms:input_files/boolean-elim10.ms
 output_files/boolean-elim10.g2.res:output_files/boolean-elim10.g2.res
 test/diff/diff_boolean-elim10.sh:test/diff/diff_boolean-elim10.sh
 input_files/boolean-random12.ms:input_files/boolean-random12.ms
 output_files/boolean-random12.g2.res:output_files/boolean-random12.g2.res
 output_files/boolean-random12.lex.g2.res:output_files/boolean-random12.lex.g2.res
 test/diff/diff_boolean-random12.sh:test/diff/diff_boolean-random12.sh
 test/diff/diff_boolean-zeros-random12.sh:test/diff/diff_boolean-zeros-random12.sh
 test/diff/diff_boolean-zeros-lex-random12.sh:test/diff/diff_boolean-zeros-lex-random12.sh
 input_files/bug-2nd-prime-bad.ms:input_files/bug-2nd-prime-bad.ms
 output_files/bug-2nd-prime-bad.g2.res:output_files/bug-2nd-prime-bad.g2.res
 test/diff/diff_bug-2nd-prime-bad.sh:test/diff/diff_bug-2nd-prime-bad.sh
 input_files/bug-68.ms:input_files/bug-68.ms
 output_files/bug-68.res:output_files/bug-68.res
 test/diff/diff_bug-68.sh:test/diff/diff_bug-68.sh
 input_files/bug-empty-tracer.ms:input_files/bug-empty-tracer.ms
 output_files/bug-empty-tracer.g2.res:output_files/bug-empty-tracer.g2.res
 test/diff/diff_bug-empty-tracer.sh:test/diff/diff_bug-empty-tracer.sh
 test/diff/diff_bug-prob-sparse-dense-la.sh:test/diff/diff_bug-prob-sparse-dense-la.sh
 input_files/choice-linear-form-qq.ms:input_files/choice-linear-form-qq.ms
 output_files/choice-linear-form-qq.P2.d0.res:output_files/choice-linear-form-qq.P2.d0.res
 output_files/choice-linear-form-qq.p128.res:output_files/choice-linear-form-qq.p128.res
 output_files/choice-linear-form-qq.p1024.res:output_files/choice-linear-form-qq.p1024.res
 test/diff/diff_choice-linear-form-qq.sh:test/diff/diff_choice-linear-form-qq.sh
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
 output_files/cyclic5-qq.P2.d4.sep0.res:output_files/cyclic5-qq.P2.d4.sep0.res
 test/diff/diff_cyclic5-qq.sh:test/diff/diff_cyclic5-qq.sh
 test/diff/diff_cyclic5-qq-options.sh:test/diff/diff_cyclic5-qq-options.sh
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
 input_files/elim-drl-31.ms:input_files/elim-drl-31.ms
 output_files/elim-drl-31.g2.e1.res:output_files/elim-drl-31.g2.e1.res
 output_files/elim-drl-31.g2.e4.res:output_files/elim-drl-31.g2.e4.res
 output_files/elim-drl-31.g2.e5.res:output_files/elim-drl-31.g2.e5.res
 test/diff/diff_elim-drl-31.sh:test/diff/diff_elim-drl-31.sh
 input_files/elim-drl-qq.ms:input_files/elim-drl-qq.ms
 output_files/elim-drl-qq.g2.e1.res:output_files/elim-drl-qq.g2.e1.res
 output_files/elim-drl-qq.g2.e4.res:output_files/elim-drl-qq.g2.e4.res
 output_files/elim-drl-qq.g2.e5.res:output_files/elim-drl-qq.g2.e5.res
 test/diff/diff_elim-drl-qq.sh:test/diff/diff_elim-drl-qq.sh
 input_files/empty-elim-qq.ms:input_files/empty-elim-qq.ms
 output_files/empty-elim-qq.g2.e1.res:output_files/empty-elim-qq.g2.e1.res
 output_files/empty-elim-qq.g2.e2.res:output_files/empty-elim-qq.g2.e2.res
 test/diff/diff_empty-elim-qq.sh:test/diff/diff_empty-elim-qq.sh
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
 input_files/fekete-valdes-n5-qq.ms:input_files/fekete-valdes-n5-qq.ms
 output_files/fekete-valdes-n5-qq.P1.d0.p128.res:output_files/fekete-valdes-n5-qq.P1.d0.p128.res
 output_files/fekete-valdes-n5-qq.P1.d0.p1024.res:output_files/fekete-valdes-n5-qq.P1.d0.p1024.res
 test/diff/diff_fekete-valdes-n5-qq.sh:test/diff/diff_fekete-valdes-n5-qq.sh
 input_files/groebner-g2.ms:input_files/groebner-g2.ms
 output_files/groebner-g2.g2.res:output_files/groebner-g2.g2.res
 test/diff/diff_groebner-g2.sh:test/diff/diff_groebner-g2.sh
 test/diff/diff_hankel-degenerate-16.sh:test/diff/diff_hankel-degenerate-16.sh
 input_files/sum-squares-system.ms:input_files/sum-squares-system.ms
 output_files/sum-squares-system.g2.res:output_files/sum-squares-system.g2.res
 test/diff/diff_sum-squares-system.sh:test/diff/diff_sum-squares-system.sh
 input_files/henrion5-qq.ms:input_files/henrion5-qq.ms
 output_files/henrion5-qq.P2.d0.res:output_files/henrion5-qq.P2.d0.res
 output_files/henrion5-qq.g2.res:output_files/henrion5-qq.g2.res
 output_files/henrion5-qq.g2.e1.res:output_files/henrion5-qq.g2.e1.res
 output_files/henrion5-qq.g2.e2.res:output_files/henrion5-qq.g2.e2.res
 test/diff/diff_henrion5-qq.sh:test/diff/diff_henrion5-qq.sh
 input_files/input-overflow-a-16.ms:input_files/input-overflow-a-16.ms
 output_files/input-overflow-a-16.P2.d0.res:output_files/input-overflow-a-16.P2.d0.res
 test/diff/diff_input-overflow-a-16.sh:test/diff/diff_input-overflow-a-16.sh
 input_files/input-overflow-b-16.ms:input_files/input-overflow-b-16.ms
 output_files/input-overflow-b-16.res:output_files/input-overflow-b-16.res
 test/diff/diff_input-overflow-b-16.sh:test/diff/diff_input-overflow-b-16.sh
 input_files/issue-230.ms:input_files/issue-230.ms
 output_files/issue-230.P1.d0.p128.res:output_files/issue-230.P1.d0.p128.res
 output_files/issue-230.P1.d0.p1024.res:output_files/issue-230.P1.d0.p1024.res
 test/diff/diff_issue-230.sh:test/diff/diff_issue-230.sh
 input_files/issue-230-squared.ms:input_files/issue-230-squared.ms
 output_files/issue-230-squared.P1.d0.p128.res:output_files/issue-230-squared.P1.d0.p128.res
 output_files/issue-230-squared.P1.d0.p1024.res:output_files/issue-230-squared.P1.d0.p1024.res
 test/diff/diff_issue-230-squared.sh:test/diff/diff_issue-230-squared.sh
 input_files/issue-310-bivariate.ms:input_files/issue-310-bivariate.ms
 output_files/issue-310-bivariate.p128.res:output_files/issue-310-bivariate.p128.res
 output_files/issue-310-bivariate.p1024.res:output_files/issue-310-bivariate.p1024.res
 test/diff/diff_issue-310-bivariate.sh:test/diff/diff_issue-310-bivariate.sh
 input_files/kat6-31.ms:input_files/kat6-31.ms
 output_files/kat6-31.P2.d0.res:output_files/kat6-31.P2.d0.res
 test/diff/diff_kat6-31.sh:test/diff/diff_kat6-31.sh
 input_files/kat6-32.ms:input_files/kat6-32.ms
 output_files/kat6-32.g2.res:output_files/kat6-32.g2.res
 test/diff/diff_kat6-32.sh:test/diff/diff_kat6-32.sh
 input_files/kat6-52.ms:input_files/kat6-52.ms
 output_files/kat6-52.g2.res:output_files/kat6-52.g2.res
 test/diff/diff_kat6-52.sh:test/diff/diff_kat6-52.sh
 input_files/kat6-gf16.ms:input_files/kat6-gf16.ms
 output_files/kat6-gf16.g2.res:output_files/kat6-gf16.g2.res
 test/diff/diff_kat6-gf16.sh:test/diff/diff_kat6-gf16.sh
 input_files/kat6-gf256.ms:input_files/kat6-gf256.ms
 output_files/kat6-gf256.g2.res:output_files/kat6-gf256.g2.res
 test/diff/diff_kat6-gf256.sh:test/diff/diff_kat6-gf256.sh
 input_files/kat6sq-qq.ms:input_files/kat6sq-qq.ms
 output_files/kat6sq-qq.x2.res:output_files/kat6sq-qq.x2.res
 test/diff/diff_kat6sq-qq.sh:test/diff/diff_kat6sq-qq.sh
 input_files/kat7-30.ms:input_files/kat7-30.ms
 output_files/kat7-30.P2.d0.res:output_files/kat7-30.P2.d0.res
 test/diff/diff_kat7-30.sh:test/diff/diff_kat7-30.sh
 input_files/kat7-qq.ms:input_files/kat7-qq.ms
 output_files/kat7-qq.P2.d0.res:output_files/kat7-qq.P2.d0.res
 test/diff/diff_kat7-qq.sh:test/diff/diff_kat7-qq.sh
 input_files/kat8-qq-truncate.ms:input_files/kat8-qq-truncate.ms
 output_files/kat8-qq-truncate.res:output_files/kat8-qq-truncate.res
 test/diff/diff_kat8-qq-truncate.sh:test/diff/diff_kat8-qq-truncate.sh
 input_files/lexgb-longchain-31.ms:input_files/lexgb-longchain-31.ms
 output_files/lexgb-longchain-31.x2.res:output_files/lexgb-longchain-31.x2.res
 test/diff/diff_lexgb-longchain-31.sh:test/diff/diff_lexgb-longchain-31.sh
 input_files/lexgb-lowerset-31.ms:input_files/lexgb-lowerset-31.ms
 output_files/lexgb-lowerset-31.x2.res:output_files/lexgb-lowerset-31.x2.res
 test/diff/diff_lexgb-lowerset-31.sh:test/diff/diff_lexgb-lowerset-31.sh
 input_files/lexgb-lowerset-qq.ms:input_files/lexgb-lowerset-qq.ms
 output_files/lexgb-lowerset-qq.x2.res:output_files/lexgb-lowerset-qq.x2.res
 test/diff/diff_lexgb-lowerset-qq.sh:test/diff/diff_lexgb-lowerset-qq.sh
 input_files/lexgb-nonradical-31.ms:input_files/lexgb-nonradical-31.ms
 output_files/lexgb-nonradical-31.x2.res:output_files/lexgb-nonradical-31.x2.res
 test/diff/diff_lexgb-nonradical-31.sh:test/diff/diff_lexgb-nonradical-31.sh
 input_files/lexgb-nonradical-qq.ms:input_files/lexgb-nonradical-qq.ms
 output_files/lexgb-nonradical-qq.x2.res:output_files/lexgb-nonradical-qq.x2.res
 test/diff/diff_lexgb-nonradical-qq.sh:test/diff/diff_lexgb-nonradical-qq.sh
 input_files/lexgb-nonshape-31.ms:input_files/lexgb-nonshape-31.ms
 output_files/lexgb-nonshape-31.x2.res:output_files/lexgb-nonshape-31.x2.res
 test/diff/diff_lexgb-nonshape-31.sh:test/diff/diff_lexgb-nonshape-31.sh
 input_files/lexgb-nonshape-qq.ms:input_files/lexgb-nonshape-qq.ms
 output_files/lexgb-nonshape-qq.x2.res:output_files/lexgb-nonshape-qq.x2.res
 test/diff/diff_lexgb-nonshape-qq.sh:test/diff/diff_lexgb-nonshape-qq.sh
 input_files/lexgb-shape-31.ms:input_files/lexgb-shape-31.ms
 output_files/lexgb-shape-31.x2.res:output_files/lexgb-shape-31.x2.res
 test/diff/diff_lexgb-shape-31.sh:test/diff/diff_lexgb-shape-31.sh
 input_files/lexgb-shape-qq.ms:input_files/lexgb-shape-qq.ms
 output_files/lexgb-shape-qq.x2.res:output_files/lexgb-shape-qq.x2.res
 test/diff/diff_lexgb-shape-qq.sh:test/diff/diff_lexgb-shape-qq.sh
 input_files/lifting-single-point.ms:input_files/lifting-single-point.ms
 output_files/lifting-single-point.p128.res:output_files/lifting-single-point.p128.res
 output_files/lifting-single-point.p1024.res:output_files/lifting-single-point.p1024.res
 test/diff/diff_lifting-single-point.sh:test/diff/diff_lifting-single-point.sh
 input_files/linear0-qq.ms:input_files/linear0-qq.ms
 output_files/linear0-qq.P2.d0.res:output_files/linear0-qq.P2.d0.res
 test/diff/diff_linear0-qq.sh:test/diff/diff_linear0-qq.sh
 input_files/linear1-qq.ms:input_files/linear1-qq.ms
 output_files/linear1-qq.P2.d0.res:output_files/linear1-qq.P2.d0.res
 test/diff/diff_linear1-qq.sh:test/diff/diff_linear1-qq.sh
 input_files/linear2-qq.ms:input_files/linear2-qq.ms
 output_files/linear2-qq.P2.d0.res:output_files/linear2-qq.P2.d0.res
 test/diff/diff_linear2-qq.sh:test/diff/diff_linear2-qq.sh
 input_files/maxbitsize-bug.ms:input_files/maxbitsize-bug.ms
 output_files/maxbitsize-bug.res:output_files/maxbitsize-bug.res
 test/diff/diff_maxbitsize-bug.sh:test/diff/diff_maxbitsize-bug.sh
 input_files/mq-2-1.ms:input_files/mq-2-1.ms
 output_files/mq-2-1.res:output_files/mq-2-1.res
 output_files/mq-2-1.P2.d0.res:output_files/mq-2-1.P2.d0.res
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
 output_files/nonradical-shape-31.P2.d0.res:output_files/nonradical-shape-31.P2.d0.res
 output_files/nonradical-shape-31.x2.res:output_files/nonradical-shape-31.x2.res
 test/diff/diff_nonradical-shape-31.sh:test/diff/diff_nonradical-shape-31.sh
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
 output_files/one-31.x1.res:output_files/one-31.x1.res
 test/diff/diff_one-31.sh:test/diff/diff_one-31.sh
 input_files/one-8.ms:input_files/one-8.ms
 output_files/one-8.res:output_files/one-8.res
 test/diff/diff_one-8.sh:test/diff/diff_one-8.sh
 input_files/one-qq.ms:input_files/one-qq.ms
 output_files/one-qq.res:output_files/one-qq.res
 output_files/one-qq.g2.res:output_files/one-qq.g2.res
 output_files/one-qq.x2.res:output_files/one-qq.x2.res
 test/diff/diff_one-qq.sh:test/diff/diff_one-qq.sh
 input_files/quadratic-nonradical-qq.ms:input_files/quadratic-nonradical-qq.ms
 output_files/quadratic-nonradical-qq.P2.d0.res:output_files/quadratic-nonradical-qq.P2.d0.res
 output_files/quadratic-nonradical-qq.p128.res:output_files/quadratic-nonradical-qq.p128.res
 output_files/quadratic-nonradical-qq.p1024.res:output_files/quadratic-nonradical-qq.p1024.res
 test/diff/diff_quadratic-nonradical-qq.sh:test/diff/diff_quadratic-nonradical-qq.sh
 input_files/radical-shape-31.ms:input_files/radical-shape-31.ms
 output_files/radical-shape-31.P2.d0.res:output_files/radical-shape-31.P2.d0.res
 output_files/radical-shape-31.x2.res:output_files/radical-shape-31.x2.res
 test/diff/diff_radical-shape-31.sh:test/diff/diff_radical-shape-31.sh
 input_files/radical-shape-qq.ms:input_files/radical-shape-qq.ms
 output_files/radical-shape-qq.P2.d0.res:output_files/radical-shape-qq.P2.d0.res
 test/diff/diff_radical-shape-qq.sh:test/diff/diff_radical-shape-qq.sh
 input_files/random6-gf256-int.ms:input_files/random6-gf256-int.ms
 output_files/random6-gf256-int.g2.res:output_files/random6-gf256-int.g2.res
 test/diff/diff_random6-gf256-int.sh:test/diff/diff_random6-gf256-int.sh
 input_files/random8-gf2.ms:input_files/random8-gf2.ms
 output_files/random8-gf2.g2.res:output_files/random8-gf2.g2.res
 test/diff/diff_random8-gf2.sh:test/diff/diff_random8-gf2.sh
 input_files/realroot1.ms:input_files/realroot1.ms
 output_files/realroot1.p0.res:output_files/realroot1.p0.res
 output_files/realroot1.p128.res:output_files/realroot1.p128.res
 output_files/realroot1.p1024.res:output_files/realroot1.p1024.res
 test/diff/diff_realroot1.sh:test/diff/diff_realroot1.sh
 input_files/realroot-extraction-exact-root.ms:input_files/realroot-extraction-exact-root.ms
 output_files/realroot-extraction-exact-root.p128.res:output_files/realroot-extraction-exact-root.p128.res
 output_files/realroot-extraction-exact-root.p1024.res:output_files/realroot-extraction-exact-root.p1024.res
 test/diff/diff_realroot-extraction-exact-root.sh:test/diff/diff_realroot-extraction-exact-root.sh
 input_files/reals-close-components.ms:input_files/reals-close-components.ms
 output_files/reals-close-components.p128.res:output_files/reals-close-components.p128.res
 output_files/reals-close-components.p1024.res:output_files/reals-close-components.p1024.res
 test/diff/diff_reals-close-components.sh:test/diff/diff_reals-close-components.sh
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
 input_files/reals-exact-left-end.ms:input_files/reals-exact-left-end.ms
 output_files/reals-exact-left-end.p128.res:output_files/reals-exact-left-end.p128.res
 output_files/reals-exact-left-end.p1024.res:output_files/reals-exact-left-end.p1024.res
 test/diff/diff_reals-exact-left-end.sh:test/diff/diff_reals-exact-left-end.sh
 input_files/reals-zero-equal.ms:input_files/reals-zero-equal.ms
 output_files/reals-zero-equal.p128.res:output_files/reals-zero-equal.p128.res
 output_files/reals-zero-equal.p1024.res:output_files/reals-zero-equal.p1024.res
 test/diff/diff_reals-zero-equal.sh:test/diff/diff_reals-zero-equal.sh
 input_files/test-lifting2.ms:input_files/test-lifting2.ms
 output_files/test-lifting2.P2.d0.res:output_files/test-lifting2.P2.d0.res
 test/diff/diff_test-lifting2.sh:test/diff/diff_test-lifting2.sh
 input_files/univariate-assign-left.ms:input_files/univariate-assign-left.ms
 output_files/univariate-assign-left.p128.res:output_files/univariate-assign-left.p128.res
 output_files/univariate-assign-left.p1024.res:output_files/univariate-assign-left.p1024.res
 test/diff/diff_univariate-assign-left.sh:test/diff/diff_univariate-assign-left.sh
 input_files/univariate-cubic-close-roots.ms:input_files/univariate-cubic-close-roots.ms
 output_files/univariate-cubic-close-roots.p128.res:output_files/univariate-cubic-close-roots.p128.res
 output_files/univariate-cubic-close-roots.p1024.res:output_files/univariate-cubic-close-roots.p1024.res
 test/diff/diff_univariate-cubic-close-roots.sh:test/diff/diff_univariate-cubic-close-roots.sh
 input_files/univariate-gcd.ms:input_files/univariate-gcd.ms
 output_files/univariate-gcd.p128.res:output_files/univariate-gcd.p128.res
 output_files/univariate-gcd.p1024.res:output_files/univariate-gcd.p1024.res
 test/diff/diff_univariate-gcd.sh:test/diff/diff_univariate-gcd.sh
 input_files/univariate-middle-root.ms:input_files/univariate-middle-root.ms
 output_files/univariate-middle-root.p128.res:output_files/univariate-middle-root.p128.res
 output_files/univariate-middle-root.p1024.res:output_files/univariate-middle-root.p1024.res
 test/diff/diff_univariate-middle-root.sh:test/diff/diff_univariate-middle-root.sh
 input_files/univariate-sgn.ms:input_files/univariate-sgn.ms
 output_files/univariate-sgn.p128.res:output_files/univariate-sgn.p128.res
 output_files/univariate-sgn.p1024.res:output_files/univariate-sgn.p1024.res
 output_files/univariate-sgn.x2.res:output_files/univariate-sgn.x2.res
 test/diff/diff_univariate-sgn.sh:test/diff/diff_univariate-sgn.sh
 input_files/univariate-wilkinson-both-signs.ms:input_files/univariate-wilkinson-both-signs.ms
 output_files/univariate-wilkinson-both-signs.p128.res:output_files/univariate-wilkinson-both-signs.p128.res
 output_files/univariate-wilkinson-both-signs.p1024.res:output_files/univariate-wilkinson-both-signs.p1024.res
 test/diff/diff_univariate-wilkinson-both-signs.sh:test/diff/diff_univariate-wilkinson-both-signs.sh
])
AC_OUTPUT
AC_MSG_NOTICE([SIMD kernels:${msolve_simd_kernels:- none}])
