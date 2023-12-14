#!/bin/bash

file=realroot-extract

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 256 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.p256.res
if [ $? -gt 0 ]; then
    exit 4
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 21
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 22
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 256 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 23
fi

diff test/diff/$file.res output_files/$file.p256.res
if [ $? -gt 0 ]; then
    exit 24
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 41
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 42
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 256 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 43
fi

diff test/diff/$file.res output_files/$file.p256.res
if [ $? -gt 0 ]; then
    exit 44
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 61
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 62
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -p 256 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 63
fi

diff test/diff/$file.res output_files/$file.p256.res
if [ $? -gt 0 ]; then
    exit 64
fi

rm test/diff/$file.res
