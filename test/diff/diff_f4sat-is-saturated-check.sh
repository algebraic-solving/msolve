#!/bin/bash

file=f4sat-is-saturated-check

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.1.res \
      -S -n 2
# should return an error 1 status
if [ $? -gt 1 ]; then
    print_exit 1
fi

diff test/diff/$file.1.res output_files/$file.res
# should return an error 1 status
if [ $? -gt 1 ]; then
    print_exit 2
fi

rm test/diff/$file.1.res

normal_exit
