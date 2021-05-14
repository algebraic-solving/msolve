lib_LTLIBRARIES 		= libneogb.la
libneogb_la_SOURCES = libneogb.h gb.c
libmsolve_ladir				=	$(includedir)/neogb
libmsolve_la_HEADERS	=	data.h f4.h
libneogb_la_LDFLAGS	= -release $(PACKAGE_VERSION)
libneogb_la_CFLAGS	= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)

EXTRA_DIST	=		basis.h \
								data.h \
								f4.h \
								hash.h \
								io.h \
								modular.h \
								nf.h \
								sort_r.h \
								stat.h \
								tools.h \
								update.h \
								basis.c \
								convert.c \
								data.c \
								f4.c \
								hash.c \
								io.c \
								la_ff_16.c \
								la_ff_32.c \
								la_ff_8.c \
								la_qq.c \
								modular.c \
								nf.c \
								order.c \
								stat.c \
								symbol.c \
								tools.c \
								update.c
