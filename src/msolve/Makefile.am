lib_LTLIBRARIES 			= libmsolve.la
libmsolve_la_SOURCES 	= msolve.c
libmsolve_la_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
libmsolve_la_LIBADD	=	../usolve/libusolve.la ../fglm/libfglm.la ../neogb/libneogb.la

bin_PROGRAMS		=	msolve
msolve_SOURCES 	= main.c
msolve_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
msolve_LDADD 		= libmsolve.la
