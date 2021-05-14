SUBDIRS = src/usolve src/fglm src/neogb src/msolve

bin_PROGRAMS		=	msolve
msolve_SOURCES 	= src/msolve/main.c
msolve_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
msolve_LDADD 		= src/neogb/libneogb.la src/fglm/libfglm.la src/usolve/libusolve.la
