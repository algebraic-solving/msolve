#!/bin/bash

file=nf-16

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -n2
if [ $? -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 4
fi

rm test/diff/$file.res
