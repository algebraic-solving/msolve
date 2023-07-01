#!/bin/bash

file=one-qq

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -g 2
if [ $? -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    exit 4
fi

rm test/diff/$file.res

