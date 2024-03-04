lib_LTLIBRARIES 			= libmsolve.la
libmsolve_la_SOURCES 	= libmsolve.c
libmsolve_ladir				=	$(includedir)/msolve/msolve
libmsolve_la_HEADERS	= msolve-data.h msolve.h
libmsolve_la_CFLAGS		= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)
libmsolve_la_LDFLAGS	= -version-info $(LT_VERSION)
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
								../crt/mpz_CRT_ui.c \
								../upolmat/nmod_mat_extra.h \
								../upolmat/nmod_mat_poly_arith.c \
								../upolmat/nmod_mat_poly_mbasis.c \
								../upolmat/nmod_mat_poly_shift.c \
								../upolmat/nmod_poly_mat_utils.c \
								../upolmat/nmod_mat_left_nullspace.c \
								../upolmat/nmod_mat_poly.h \
								../upolmat/nmod_mat_poly_mem.c \
								../upolmat/nmod_poly_mat_pmbasis.c \
								../upolmat/nmod_poly_mat_utils.h \
								../upolmat/nmod_mat_permute_rows.c \
								../upolmat/nmod_mat_poly_io.c \
								../upolmat/nmod_mat_poly_set_from.c \
								../upolmat/nmod_poly_mat_pmbasis.h
