AUTOMAKE_OPTIONS = color-tests
SUBDIRS = src/usolve src/fglm src/neogb src/msolve

AM_CFLAGS	= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
LDADD 		= src/neogb/libneogb.la src/fglm/libfglm.la src/usolve/libusolve.la

bin_PROGRAMS	= msolve
msolve_SOURCES 	= src/msolve/main.c

check_PROGRAMS		= neogb_io \
			  fglm_build_matrixn_radical-shape-31 \
			  fglm_build_matrixn_nonradical-shape-31 \
			  fglm_build_matrixn_nonradical-radicalshape-31 \
			  line_endings_support

checkdiff               = test/diff/diff_bug-2nd-prime-bad.sh \
			  test/diff/diff_bug-68.sh \
			  test/diff/diff_bug-empty-tracer.sh \
			  test/diff/diff_cp-d3-n4-p2.sh \
			  test/diff/diff_cyclic5-16.sh \
			  test/diff/diff_cyclic5-31.sh \
			  test/diff/diff_cyclic5-qq.sh \
			  test/diff/diff_eco10-31.sh \
			  test/diff/diff_eco6-16.sh \
			  test/diff/diff_eco6-31.sh \
			  test/diff/diff_eco6-qq.sh \
			  test/diff/diff_elim-31.sh \
			  test/diff/diff_elim-qq.sh \
			  test/diff/diff_empty-elim-qq.sh \
			  test/diff/diff_f4sat-31.sh \
			  test/diff/diff_f4sat-byone-31.sh \
			  test/diff/diff_f4sat-field-char.sh \
			  test/diff/diff_f4sat-is-saturated-check.sh \
			  test/diff/diff_f4sat-zero-input.sh \
			  test/diff/diff_groebner-g2.sh \
			  test/diff/diff_henrion5-qq.sh \
			  test/diff/diff_input-overflow-a-16.sh \
			  test/diff/diff_input-overflow-b-16.sh \
			  test/diff/diff_issue-230.sh \
			  test/diff/diff_issue-230-squared.sh \
			  test/diff/diff_kat6-31.sh \
			  test/diff/diff_kat7-qq.sh \
			  test/diff/diff_kat8-qq-truncate.sh \
			  test/diff/diff_linear0-qq.sh \
			  test/diff/diff_linear1-qq.sh \
			  test/diff/diff_linear2-qq.sh \
			  test/diff/diff_maxbitsize-bug.sh \
			  test/diff/diff_mq-2-1.sh \
			  test/diff/diff_multy-16.sh \
			  test/diff/diff_multy-31.sh \
			  test/diff/diff_multy-qq.sh \
			  test/diff/diff_nf-16.sh \
			  test/diff/diff_nf-31.sh \
			  test/diff/diff_nf-8.sh \
			  test/diff/diff_nf-lm-bug.sh \
			  test/diff/diff_nonradical-radicalshape-no-square-31.sh \
			  test/diff/diff_nonradical-radicalshape-no-square-qq.sh \
			  test/diff/diff_nonradical-radicalshape-qq.sh \
			  test/diff/diff_nonradical-shape-31.sh \
			  test/diff/diff_nonradical-shape-qq.sh \
			  test/diff/diff_one-16.sh \
			  test/diff/diff_one-31.sh \
			  test/diff/diff_one-qq.sh \
			  test/diff/diff_quadratic-nonradical-qq.sh \
			  test/diff/diff_radical-shape-31.sh \
			  test/diff/diff_radical-shape-qq.sh \
			  test/diff/diff_realroot1.sh \
			  test/diff/diff_realroot-extraction-exact-root.sh \
			  test/diff/diff_reals-dim0-chgvar.sh \
			  test/diff/diff_reals-dim0-extract.sh \
			  test/diff/diff_reals-dim0.sh \
			  test/diff/diff_reals-dim0-swapvar.sh \
			  test/diff/diff_test-lifting2.sh

# dist_check_DATA         = test/input_files
neogb_io_SOURCES 	= test/neogb/io/validate_input_data.c
fglm_build_matrixn_radical_shape_31_SOURCES = test/fglm/build_matrixn_radical-shape-31.c
fglm_build_matrixn_nonradical_shape_31_SOURCES = test/fglm/build_matrixn_nonradical-shape-31.c
fglm_build_matrixn_nonradical_radicalshape_31_SOURCES = test/fglm/build_matrixn_nonradical-radicalshape-31.c
line_endings_support_SOURCES = test/msolve/line_endings_support.c

TESTS = $(check_PROGRAMS) $(checkdiff)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = msolve.pc
