#!/bin/bash

file=nonradical_shape-qq

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -p 1024
if [ ! $? ]; then
    exit 3
fi

diff test/diff/$file.p1024.res output_files/$file.res
if [ !  $? ]; then
    exit 4
fi

rm test/diff/$file.res
