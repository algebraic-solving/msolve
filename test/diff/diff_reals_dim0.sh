#!/bin/bash

file=reals_dim0

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res
if [ ! $? ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ !  $? ]; then
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -p 256
if [ ! $? ]; then
    exit 3
fi

diff test/diff/$file.p256.res output_files/$file.res
if [ !  $? ]; then
    exit 4
fi

rm test/diff/$file.res
