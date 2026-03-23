#!/usr/bin/env bash

file=f4sat-zero-input

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -S -g 2 -n 2
# should return an error 1 status
if [ $? -gt 1 ]; then
    print_exit 1
fi

diff test/diff/$file.res output_files/$file.res
# should return an error 1 status
if [ $? -gt 1 ]; then
    print_exit 2
fi

rm test/diff/$file.res

normal_exit
