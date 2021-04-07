CC = gcc 

$(eval GIT_COMMIT_SHA ="$(shell git rev-parse HEAD)")
GIT_C_SHA := -D'GIT_COMMIT_HASH=$(GIT_COMMIT_SHA)'

LOCAL_CFLAGS = -std=c99 -O3 -DNDEBUG -DUSEFLINT -Wall -Winline  -funroll-loops -fomit-frame-pointer -Wall -Wshadow -Wpointer-arith -Wwrite-strings -static-libgcc -fPIC #$(LIBGOMP)

LIBGOMP = /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.a
LIBS = -I/usr/local/include -L/usr/local/lib/ -I$(CURDIR) 

RM = rm

GMP_INC_DIR=.
GMP_LIB_DIR=.

LIB_PATHS = -L/usr/local/lib
INC_PATHS = -I/usr/local/include

# -L/usr/local/lib/libgmp.a -L/usr/local/lib/libflint.a

CFLAGS = -c

LDFLAGS_F = -L/usr/local/lib/libmpfr.a -L/usr/local/lib/libgmp.a -L/usr/local/lib/libflint.a -fopenmp -lpthread -lflint -lmpfr -lgmp  -lm
LDFLAGS = -Wall $(DEBUG) $(LIB_PATHS) $(INC_PATHS) -lmpfr -lgmp -lpthread -fopenmp -lm


all : tomaple tomathematica usolve

usolve : main.c usolve.c utils.c taylor_shift.c descartes.c evaluate.c print_usolve.c debug.c refine.c
	 $(CC) -static -D'GIT_COMMIT_HASH=$(GIT_COMMIT_SHA)' $(LOCAL_CFLAGS) -o usolve main.c $(LIB_PATHS) $(INC_PATHS) $(LDFLAGS_F)

libusolve.o : usolve.c utils.c taylor_shift.c descartes.c evaluate.c print_usolve.c debug.c refine.c
	$(CC) -c -static usolve.c $(LOCAL_CFLAGS) $(LDFLAGS_F) $(LIB_PATHS) $(INC_PATHS) -o libusolve.o
		ar rcsv libusolve.a /usr/local/lib/libmpfr.a /usr/local/lib/libgmp.a /usr/local/lib/libflint.a libusolve.o
		ranlib libusolve.a

tomaple : tomaple.o
	$(CC) -static -o tomaple tomaple.o   $(LDFLAGS)

tomaple.o : tomaple.c
	$(CC) -c tomaple.c $(LDFLAGS)


tomathematica : tomathematica.o
	$(CC) -static -o tomathematica tomathematica.o  /usr/local/lib/libmpfr.a /usr/local/lib/libgmp.a  $(LDFLAGS)


tomathematica.o : tomathematica.c
	$(CC) -c tomathematica.c $(LDFLAGS)

# taylor_test : utils.o univmultiply.o mpz_upoly_multiply.o taylor_shift.o tests.c
# 	$(CC) -static $(LOCAL_CFLAGS) -o taylor_test tests.c utils.o univmultiply.o mpz_upoly_multiply.o taylor_shift.o $(LDFLAGS) 


clean:
	rm usolve libusolve.o

.PHONY: all clean
