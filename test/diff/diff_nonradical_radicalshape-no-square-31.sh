#!/bin/bash

file=nonradical_radicalshape-no-square-31

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 2 -t 1 -c 0
if [ $?  -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    exit 4
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 2 -t 2 -c 0
if [ $?  -gt 0 ]; then
    exit 23
fi

diff test/diff/$file.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    exit 24
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 44 -t 1 -c 0
if [ $?  -gt 0 ]; then
    exit 43
fi

diff test/diff/$file.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    exit 44
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 44 -t 2 -c 0
if [ $?  -gt 0 ]; then
    exit 63
fi

diff test/diff/$file.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    exit 64
fi

rm test/diff/$file.res
