AUTOMAKE_OPTIONS = color-tests
SUBDIRS = src/usolve src/fglm src/neogb src/msolve/

AM_CFLAGS	= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
LDADD 		= src/neogb/libneogb.la src/fglm/libfglm.la src/usolve/libusolve.la

bin_PROGRAMS	= msolve
msolve_SOURCES 	= src/msolve/main.c

check_PROGRAMS		= neogb_io fglm_build_matrixn
neogb_io_SOURCES 	= test/neogb/io/validate_input_data.c
fglm_build_matrixn_SOURCES            = test/fglm/build_matrixn.c

TESTS = $(check_PROGRAMS)
