#!/bin/bash

file=nf-8

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -n2
if [ $? -gt 0 ]; then
    exit 5
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 6
fi

rm test/diff/$file.res
