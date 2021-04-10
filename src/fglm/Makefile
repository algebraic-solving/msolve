# gcc -Wall -O3 -fopenmp -fprofile-generate -fopt-info-vec-all -ffast-math -funroll-all-loops -ftree-vectorize  -mavx2 -march=native main.c 2> /tmp/compile
# /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.a

CC = gcc -O2 -Wall -funroll-loops -fomit-frame-pointer -fexpensive-optimizations -march=native -fPIC

LIB_PATHS = -L/usr/local/lib/
INC_PATHS = -I/usr/local/include 

CFLAGS = -Wall -Wshadow -Wpointer-arith -Wwrite-strings

LDFLAGS = -Wall $(DEBUG) $(LIB_PATHS) -static-libgcc -fopenmp -lpthread -lflint -lmpfr -lgmp -lm
LDFLAGS0 = -Wall $(DEBUG) $(LIB_PATHS) -lgmp  -lmpfr -lflint  -lm 

fglm: main.c linalg-fglm.c berlekamp_massey.c data_fglm.c fglm_core.c fglm_core.o
		$(CC) -static main.c -o fglm $(LDFLAGS) /usr/local/lib/libmpfr.a /usr/local/lib/libgmp.a /usr/local/lib/libflint.a #/usr/lib/gcc/x86_64-linux-gnu/5/libgomp.a

libfglm.o: linalg-fglm.c berlekamp_massey.c data_fglm.c fglm_core.c
		$(CC) -c -static fglm_core.c $(LDFLAGS)  -o libfglm.o

libfglm.a: libfglm.o
		ar rcsv libfglm.a /usr/local/lib/libmpfr.a /usr/local/lib/libgmp.a /usr/local/lib/libflint.a  libfglm.o
		ranlib libfglm.a

all:

clean:
	rm libfglm.o 
