#!/bin/bash

file=nonradical_radicalshape-no-square-31

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.3.res \
      -P 2 -l 2 -t 1 -c 0
if [ $?  -gt 0 ]; then
    print_exit 3
fi

diff test/diff/$file.3.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    print_exit 4
fi

rm test/diff/$file.3.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.23.res \
      -P 2 -l 2 -t 2 -c 0
if [ $?  -gt 0 ]; then
    print_exit 23
fi

diff test/diff/$file.23.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    print_exit 24
fi

rm test/diff/$file.23.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.43.res \
      -P 2 -l 44 -t 1 -c 0
if [ $?  -gt 0 ]; then
    print_exit 43
fi

diff test/diff/$file.43.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    print_exit 44
fi

rm test/diff/$file.43.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.63.res \
      -P 2 -l 44 -t 2 -c 0
if [ $?  -gt 0 ]; then
    print_exit 63
fi

diff test/diff/$file.63.res output_files/$file.c0.res
if [ $? -gt 0 ]; then
    print_exit 64
fi

rm test/diff/$file.63.res
