CC = gcc
local_CFLAGS = -std=c99 -O2 -Wall -pedantic -funroll-loops -fomit-frame-pointer -fexpensive-optimizations -march=native -fPIC -static-libgcc $(CFLAGS) $(LIBGOMP)

LIBGOMP = /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.a
#LIBGOMP = /usr/local/lib64/libgomp.a

LIBS = -I/usr/local/include -L/usr/local/lib/ -I$(CURDIR) 

ALLSOURCES := $(wildcard *.c)

LDFLAGS = $(LIBS) -lmpfr -lgmp -lpthread -fopenmp -lm

$(eval GIT_COMMIT_SHA ="$(shell git rev-parse HEAD)")

all: msolve

.c.o:
		@echo [Compile] $<
		@$(CC) -c $(local_CFLAGS) $< -o $@ 

msolve: $(ALLSOURCES) 
		$(MAKE) libusolve.o -C ../usolve
		$(MAKE) libfglm.o -C ../fglm
		$(MAKE) libneogb.o -C ../neogb
		$(CC) -D'GIT_COMMIT_HASH=$(GIT_COMMIT_SHA)' $(local_CFLAGS) -static main.c -o msolve $(LDFLAGS) -I/usr/local/include/ -I../usolve/ ../usolve/libusolve.o -I../fglm/ ../fglm/libfglm.o -I../neogb/ ../neogb/libneogb.o -L/usr/local/lib/libmpfr.a -L$(LIBGOMP) -L/usr/local/lib/libgmp.a -L/usr/local/lib/libflint.a -lpthread -fopenmp -lgomp -lflint -lmpfr -lgmp -lm

libmsolve.so: $(ALLSOURCES)
		$(MAKE) libusolve.o -C ../usolve
		$(MAKE) libfglm.o -C ../fglm
		$(MAKE) libneogb.o -C ../neogb
		$(CC) -D'GIT_COMMIT_HASH=$(GIT_COMMIT_SHA)' $(local_CFLAGS) -shared -Wl,-soname,libmsolve.so.0 $(LDFLAGS) -o libmsolve.so.0.2.0 ../usolve/libusolve.o ../fglm/libfglm.o ../neogb/libneogb.o libmsolve.c

clean:
		rm -f *.o	
		rm -f ../usolve/*.o	
		rm -f ../fglm/*.o	
		rm -f ../neogb/*.o	
		rm -f libmsolve.a
		rm -f libmsolve.so*
		rm -f msolve

.PHONY: all check clean
