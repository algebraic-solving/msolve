MPFRDIR = /usr/local/lib/libmpfr.a
GMPDIR = /usr/local/lib/libgmp.a
FLINTDIR = /usr/local/lib/libflint.a

CC = gcc

CFLAGS = -O3 -Wall -funroll-loops -Wshadow -Wpointer-arith -Wwrite-strings -Winline -static-libgcc -fopenmp

LDFLAGS = -L$(MPFRDIR) -L$(GMPDIR) -L$(FLINTDIR)
LDLIBS = -lpthread -lflint -lmpfr -lgmp  -lm

all : libusolve.o

libusolve.o : usolve.c utils.c taylor_shift.c descartes.c evaluate.c print_usolve.c debug.c refine.c univmultiply.c data_usolve.c
	$(CC) -c -static -DNDEBUG -DUSEFLINT usolve.c $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o libusolve.o

clean:
	rm libusolve.o

.PHONY: all clean
