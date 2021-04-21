MPFRDIR = /usr/local/lib/libmpfr.a
GMPDIR = /usr/local/lib/libgmp.a
FLINTDIR = /usr/local/lib/libflint.a
LIBGOMP = /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.a
#LIBGOMP = /usr/local/lib64/libgomp.a

CC = gcc
LDLIBS = -lpthread -lgomp -lflint -lmpfr -lgmp  -lm

local_CFLAGS = -std=c99 -O2 -Wall -pedantic -funroll-loops -fomit-frame-pointer -march=native -fopenmp -fPIC -static-libgcc $(CFLAGS)

ALLSOURCES := $(wildcard *.c)

$(eval GIT_COMMIT_SHA ="$(shell git rev-parse HEAD)")

all: msolve libmsolve.so

.c.o:
		@echo [Compile] $<
		@$(CC) -c $(local_CFLAGS) $< -o $@ 

msolve: $(ALLSOURCES) 
		$(MAKE) libusolve.o -C ../usolve
		$(MAKE) libfglm.o -C ../fglm
		$(MAKE) libneogb.o -C ../neogb
		$(CC) -D'GIT_COMMIT_HASH=$(GIT_COMMIT_SHA)' $(local_CFLAGS) -static main.c -o msolve -I/usr/local/include/ -I../usolve/ ../usolve/libusolve.o -I../fglm/ ../fglm/libfglm.o -I../neogb/ ../neogb/libneogb.o $(LDLIBS) 

deps:
		$(MAKE) libusolve.o -C ../usolve
		$(MAKE) libfglm.o -C ../fglm
		$(MAKE) libneogb.o -C ../neogb

libmsolve.so: deps $(ALLSOURCES)
		$(CC) -D'GIT_COMMIT_HASH=$(GIT_COMMIT_SHA)' $(local_CFLAGS) -shared -Wl,-soname,libmsolve.so.1 $(LDFLAGS) -o libmsolve.so.0.1.1 ../usolve/libusolve.o ../fglm/libfglm.o ../neogb/libneogb.o libmsolve.c

clean:
		rm -f *.o
		rm -f ../usolve/*.o
		rm -f ../fglm/*.o
		rm -f ../neogb/*.o
		rm -f libmsolve.a
		rm -f libmsolve.so*
		rm -f msolve

.PHONY: all check clean deps
