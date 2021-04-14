CC = gcc

CFLAGS = -O3 -Wall -funroll-loops -march=native -Wshadow -Wpointer-arith -Wwrite-strings -Winline -static-libgcc -fopenmp

LDLIBS = -lpthread

all : libfglm.o

libfglm.o: linalg-fglm.c berlekamp_massey.c data_fglm.c fglm_core.c
		$(CC) -c -static fglm_core.c $(CFLAGS) $(LDLIBS)  -o libfglm.o

clean:
	rm libfglm.o
