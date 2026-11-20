#!/usr/bin/env bash

#test file for real root extraction with adaptative precision
#For this test file, one coordinate of one solution is very close to 0
#msolve's real root extraction should be able to make the distinction

#test file for the real roots with the precision -p 0: the coordinates
#are computed to 2^0.

file=realroot1

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.15.res \
      --random-seed $seed \
      -p 0 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    print_exit 15
fi

diff test/diff/$file.15.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 16
fi

rm test/diff/$file.15.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.35.res \
      --random-seed $seed \
      -p 0 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    print_exit 35
fi

diff test/diff/$file.35.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 36
fi

rm test/diff/$file.35.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.55.res \
      --random-seed $seed \
      -p 0 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    print_exit 55
fi

diff test/diff/$file.55.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 56
fi

rm test/diff/$file.55.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.75.res \
      --random-seed $seed \
      -p 0 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    print_exit 75
fi

diff test/diff/$file.75.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 76
fi

rm test/diff/$file.75.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.95.res \
      --random-seed $seed \
      -p 0 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    print_exit 95
fi

diff test/diff/$file.95.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 96
fi

rm test/diff/$file.95.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.115.res \
      --random-seed $seed \
      -p 0 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    print_exit 115
fi

diff test/diff/$file.115.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 116
fi

rm test/diff/$file.115.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.135.res \
      --random-seed $seed \
      -p 0 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    print_exit 135
fi

diff test/diff/$file.135.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 136
fi

rm test/diff/$file.135.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.155.res \
      --random-seed $seed \
      -p 0 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    print_exit 155
fi

diff test/diff/$file.155.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 156
fi

rm test/diff/$file.155.res

normal_exit
