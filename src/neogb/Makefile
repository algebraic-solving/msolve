CC = gcc
local_CFLAGS = -std=c99 -O2 -Wall -pedantic -funroll-loops -fomit-frame-pointer -fexpensive-optimizations -march=native -fPIC -static-libgcc $(CFLAGS)

LIBS = -I/usr/local/include -L/usr/local/lib/ -I$(CURDIR) 
LIBGOMP = /usr/lib/gcc/x86_64-linux-gnu/9/
#LIBGOMP = /usr/local/lib64/

UNITTESTS	:= $(wildcard $(CURDIR)/../../tests/neogb/*.c)
TRGUTESTS := $(UNITTESTS:.c=)
BENCHTESTS:= $(wildcard $(CURDIR)/../../benchs/neogb/*.c)
TRGBTESTS := $(BENCHTESTS:.c=)
UPROGSTST := $(patsubst %.c,%,$(UNITTESTS))
UNITPROGS	:= $(sort $(UPROGSTST))
BPROGSTST := $(patsubst %.c,%,$(BENCHTESTS))
BENCHPROGS:= $(sort $(BPROGSTST))
ALLSOURCES := $(wildcard *.c)

STATIC	=	libneogb.a
SRC	=	gb.c

OBJ = $(SRC:.c=.o)

LDFLAGS = $(LIBS) -lmpfr -lgmp -lpthread -fopenmp -lm

all: libneogb.o neogb

libneogb.o: libneogb.h gb.c
		$(CC) $(local_CFLAGS) -c -static gb.c $(LDFLAGS) -o libneogb.o
		ar rcsv libneogb.a /usr/local/lib/libmpfr.a /usr/local/lib/libgmp.a $(LIBGOMP)/libgomp.a libneogb.o
		ranlib libneogb.a

# $(STATIC): $(OBJ) data.h
#     @echo "[Link (Static)]"
#     @ar rcv $@ $^
#     @ranlib $(STATIC)

.c.o:
		@echo [Compile] $<
		@$(CC) -c $(local_CFLAGS) $< -o $@

neogb: $(ALLSOURCES) 
		$(CC) $(local_CFLAGS) -static main.c -o neogb $(LDFLAGS) -I/usr/local/include/ -I../fglm/ ../fglm/libfglm.o -L/usr/local/lib/libmpfr.a -L/usr/local/lib64/libgomp.a -L/usr/local/lib/libgmp.a -L/usr/local/lib/libmpir.a -L/usr/local/lib/libflint.a -lpthread -fopenmp -lgomp -lflint -lmpfr -lgmp -lm

check: 
		for f in $(UNITPROGS) ; do \
				$(CC) $(local_CFLAGS) -o $$f $$f.c $(LDFLAGS) $(CURDIR)/$(STATIC); \
		done
		for f in $(UNITPROGS) ; do \
				$(CURDIR)/test-driver --color-tests yes \
				--test-name "$$f" --log-file $$f.log \
				--trs-file $$f.trs $$f; \
		done

bench:
		for f in $(BENCHPROGS) ; do \
				$(CC) $(local_CFLAGS) -o $$f $$f.c $(LDFLAGS) $(CURDIR)/$(STATIC); \
		done

clean:
		rm -f *.o	
		rm -f libneogb.a
		rm -f neogb
		rm -f *.pbm
		rm -f $(TRGUTESTS)
		rm -f $(CURDIR)/../../tests/neogb/*.log
		rm -f $(CURDIR)/../../tests/neogb/*.trs
		rm -f $(CURDIR)/../../tests/neogb/*.pbm
		rm -f $(TRGBTESTS)

.PHONY: all check clean
