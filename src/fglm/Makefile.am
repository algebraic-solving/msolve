noinst_LTLIBRARIES 	= libfglm.la
libfglm_la_SOURCES 	= fglm_core.c
libfglm_ladir				=	$(includedir)/msolve/fglm
libfglm_la_HEADERS	= fglm.h
libfglm_la_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)

EXTRA_DIST	=		fglm.h \
                libfglm.h \
								berlekamp_massey.c \
								data_fglm.c \
								fglm_core.c \
								inner-product.c \
								matrix-mult.c \
								linalg-fglm.c
