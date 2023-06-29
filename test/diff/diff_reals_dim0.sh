#!/bin/bash

file=reals_dim0

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -p 256
if [ $? -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.p256.res
if [ $? -gt 0 ]; then
    exit 4
fi

rm test/diff/$file.res
