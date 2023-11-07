#!/bin/bash

file=nonradical_radicalshape-no-square-qq

# source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 3 -l 2 -t 1
if [ $?  -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.p3.res
if [ $? -gt 0 ]; then
    exit 4
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 2 -t 1 -c 0
if [ $? -ne 1 ]; then
    exit 5
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 3 -l 2 -t 2
if [ $?  -gt 0 ]; then
    exit 23
fi

diff test/diff/$file.res output_files/$file.p3.res
if [ $? -gt 0 ]; then
    exit 24
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 2 -t 2 -c 0
if [ $?  -ne 1 ]; then
    exit 25
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 3 -l 44 -t 1
if [ $?  -gt 0 ]; then
    exit 43
fi

diff test/diff/$file.res output_files/$file.p3.res
if [ $? -gt 0 ]; then
    exit 44
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 44 -t 1 -c 0
if [ $?  -ne 1 ]; then
    exit 45
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 3 -l 44 -t 2
if [ $?  -gt 0 ]; then
    exit 63
fi

diff test/diff/$file.res output_files/$file.p3.res
if [ $? -gt 0 ]; then
    exit 64
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -l 44 -t 2 -c 0
if [ $?  -ne 1 ]; then
    exit 65
fi

rm test/diff/$file.res
