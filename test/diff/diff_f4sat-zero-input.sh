#!/bin/bash

file=f4sat-zero-input

$(pwd)/msolve -f input_files/$file.ms -S -g2 -o test/diff/$file.res \
      -n2
# should return an error 1 status
if [ $? -gt 1 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
# should return an error 1 status
if [ $? -gt 1 ]; then
    exit 2
fi

rm test/diff/$file.res
