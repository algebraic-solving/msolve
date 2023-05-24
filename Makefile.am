AUTOMAKE_OPTIONS = color-tests
SUBDIRS = src/usolve src/fglm src/neogb src/msolve/
EXTRA_DIST         = $(srcdir)/test/input_files

AM_CFLAGS	= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
LDADD 		= src/neogb/libneogb.la src/fglm/libfglm.la src/usolve/libusolve.la

bin_PROGRAMS	= msolve
msolve_SOURCES 	= src/msolve/main.c

check_PROGRAMS		= neogb_io \
			  fglm_build_matrixn_radical_shape-31 \
			  fglm_build_matrixn_nonradical_shape-31 \
			  fglm_build_matrixn_nonradical_radicalshape-31
# dist_check_DATA         = test/input_files
neogb_io_SOURCES 	= test/neogb/io/validate_input_data.c
fglm_build_matrixn_radical_shape_31_SOURCES = test/fglm/build_matrixn_radical_shape-31.c
fglm_build_matrixn_nonradical_shape_31_SOURCES = test/fglm/build_matrixn_nonradical_shape-31.c
fglm_build_matrixn_nonradical_radicalshape_31_SOURCES = test/fglm/build_matrixn_nonradical_radicalshape-31.c


TESTS = $(check_PROGRAMS)
