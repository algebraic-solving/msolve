#!/bin/bash

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -d 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -L 1 -P 2 -d 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 101
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 201
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -d 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 21
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 22
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -L 1 -P 2 -d 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 211
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 221
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -d 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 41
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 42
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -L 1 -P 2 -d 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 411
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 421
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 2 -d 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 61
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 62
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -L 1 -P 2 -d 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 611
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 621
fi

rm test/diff/$file.res
