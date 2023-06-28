#!/bin/bash

file=elim-31

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -e 1 -g 2
if [ ! $? ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ !  $? ]; then
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -e 2 -g 2
if [ ! $? ]; then
    exit 3
fi

diff test/diff/$file.e2.res output_files/$file.res
if [ !  $? ]; then
    exit 4
fi

rm test/diff/$file.res
