#!/bin/bash

file=f4sat-byone-31

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.1.res \
      -S -g 2 -l 2 -t 1
if [ $? -gt 0 ]; then
    print_exit 1
fi

diff test/diff/$file.1.res output_files/$file.res
if [ $? -gt 0 ]; then
    print_exit 2
fi

rm test/diff/$file.1.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.21.res \
      -S -g 2 -l 2 -t 2
if [ $? -gt 0 ]; then
    print_exit 21
fi

diff test/diff/$file.21.res output_files/$file.res
if [ $? -gt 0 ]; then
    print_exit 22
fi

rm test/diff/$file.21.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.41.res \
      -S -g 2 -l 44 -t 1
if [ $? -gt 0 ]; then
    print_exit 41
fi

diff test/diff/$file.41.res output_files/$file.res
if [ $? -gt 0 ]; then
    print_exit 42
fi

rm test/diff/$file.41.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.61.res \
      -S -g 2 -l 44 -t 2
if [ $? -gt 0 ]; then
    print_exit 61
fi

diff test/diff/$file.61.res output_files/$file.res
if [ $? -gt 0 ]; then
    print_exit 62
fi

rm test/diff/$file.61.res

normal_exit
