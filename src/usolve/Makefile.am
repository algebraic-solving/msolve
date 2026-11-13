noinst_LTLIBRARIES 		= libusolve.la
libusolve_la_SOURCES 	= usolve.c
libusolve_ladir				=	$(includedir)/msolve/usolve
libusolve_la_HEADERS	= libusolve.h
libusolve_la_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS) -Wall -Wextra

EXTRA_DIST = 		libusolve.h \
								bisection_balanced.c \
								data_usolve.c \
								descartes.c \
								evaluate.c \
								refine.c \
								refine_arb.c \
								taylor_shift.c \
								univmultiply.c \
								utils.c
