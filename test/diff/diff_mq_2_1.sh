#!/bin/bash

file=mq_2_1

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 2 -t 1
if [ $? -gt 0 ]; then
    echo "1"
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    echo "2"
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 2 -t 2
if [ $? -gt 0 ]; then
    echo "21"
    exit 21
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    echo "22"
    exit 22
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 44 -t 1
if [ $? -gt 0 ]; then
    echo "41"
    exit 41
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    echo "42"
    exit 42
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -l 44 -t 2
if [ $? -gt 0 ]; then
    echo "61"
    exit 61
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    echo "62"
    exit 62
fi

rm test/diff/$file.res
