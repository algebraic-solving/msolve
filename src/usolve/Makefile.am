noinst_LTLIBRARIES 		= libusolve.la
libusolve_la_SOURCES 	= usolve.c
libusolve_ladir				=	$(includedir)/msolve/usolve
libusolve_la_HEADERS	= libusolve.h
libusolve_la_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS) -Wall -Wextra

EXTRA_DIST = 		libusolve.h \
                bisection.h \
								debug.h \
								descartes.h \
								evaluate.h \
								libusolve.h \
								print_usolve.h \
								refine.h \
								taylor_shift.h \
								univmultiply.h \
								bisection_balanced.c \
								data_usolve.c \
								debug.c \
								descartes.c \
								evaluate.c \
								print_usolve.c \
								refine.c \
								refine_arb.c \
								taylor_shift.c \
								tests.c \
								univmultiply.c \
								utils.c
