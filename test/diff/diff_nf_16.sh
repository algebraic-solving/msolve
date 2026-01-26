#!/bin/bash

file=nf-16

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.1.res \
      -n 2
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.1.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 2
fi

rm test/diff/$file.1.res
