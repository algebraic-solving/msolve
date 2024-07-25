#!/bin/bash

file=f4sat-field-char

$(pwd)/msolve -f input_files/$file.ms -S -g2 -o test/diff/$file.res \
      -n2
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 2
fi

rm test/diff/$file.res
