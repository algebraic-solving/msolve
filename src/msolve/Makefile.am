lib_LTLIBRARIES 			= libmsolve.la
libmsolve_la_SOURCES 	= libmsolve.c
libmsolve_ladir				=	$(includedir)/msolve/msolve
libmsolve_la_HEADERS	=
libmsolve_la_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
libmsolve_la_LDFLAGS	= -release $(PACKAGE_VERSION)
libmsolve_la_LIBADD		=	../usolve/libusolve.la ../fglm/libfglm.la ../neogb/libneogb.la

EXTRA_DIST	=	  msolve-data.h \
								msolve.h \
								msolve-data.c \
								duplicate.c \
								hilbert.c \
								linear.c \
								lifting.c \
								lifting-gb.c \
								iofiles.c \
								msolve.c \
								primes.c \
								../crt/longlong.h \
								../crt/ulong_extras.h \
								../crt/mpq_reconstruct.c \
								../crt/mpz_CRT_ui.c
