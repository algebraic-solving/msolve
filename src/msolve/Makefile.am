lib_LTLIBRARIES = libmsolve.la
libmsolve_la_SOURCES = linalg-fglm.c berlekamp_massey.c data_fglm.c fglm_core.c
libmsolve_la_LIBADD	=	../usolve/libusolve.la ../fglm/libfglm.la ../neogb/libneogb.la
bin_PROGRAMS	=	msolva
msolva_SOURCES = main.c
msolva_LDADD = libmsolve.la
