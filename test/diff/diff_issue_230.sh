#!/bin/bash

#test file for real root extraction with adaptative precision
#For this test file, one coordinate of one solutions is very close to 0
#msolve's real root extraction should be able to make the distinction

file=issue_230

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 21
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 22
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 41
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 42
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -P 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 61
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 62
fi

rm test/diff/$file.res
