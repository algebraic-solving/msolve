#!/bin/bash

file=nonradical_radicalshape-qq

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -p 3
if [ $?  -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.p3.res
if [ $? -gt 0 ]; then
    exit 4
fi

rm test/diff/$file.res
