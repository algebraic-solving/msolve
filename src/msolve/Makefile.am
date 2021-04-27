lib_LTLIBRARIES 			= libmsolve.la
libmsolve_la_SOURCES 	= libmsolve.c
libmsolve_la_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
libmsolve_la_LDFLAGS	= -release $(PACKAGE_VERSION)
libmsolve_la_LIBADD		=	../usolve/libusolve.la ../fglm/libfglm.la ../neogb/libneogb.la
