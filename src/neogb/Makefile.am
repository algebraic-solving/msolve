lib_LTLIBRARIES 		= libneogb.la
libneogb_la_SOURCES = libneogb.h gb.c
libneogb_ladir			=	$(includedir)/msolve/neogb
libneogb_la_HEADERS	=libneogb.h basis.h data.h engine.h f4.h sba.h hash.h \
					 io.h modular.h nf.h f4sat.h sort_r.h meta_data.h \
					 tools.h update.h
libneogb_la_LDFLAGS	= -release $(PACKAGE_VERSION)
libneogb_la_CFLAGS	= $(SIMD_FLAGS) $(CPUEXT_FLAGS) $(OPENMP_CFLAGS)

EXTRA_DIST	=		basis.h \
								data.h \
								engine.h \
								f4.h \
								sba.h \
								hash.h \
								io.h \
								modular.h \
								nf.h \
								f4sat.h \
								sort_r.h \
								meta_data.h \
								tools.h \
								update.h \
								basis.c \
								convert.c \
								data.c \
								engine.c \
								f4.c \
								f4sat.c \
								sba.c \
								hash.c \
								io.c \
								la_ff_16.c \
								la_ff_32.c \
								la_ff_8.c \
								la_qq.c \
								modular.c \
								nf.c \
								order.c \
								meta_data.c \
								symbol.c \
								tools.c \
								update.c
