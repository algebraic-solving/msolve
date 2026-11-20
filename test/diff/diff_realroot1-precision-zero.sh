#!/usr/bin/env bash

#test file for the real roots with the precision -p 0: the coordinates
#are computed to 2^0.

file=realroot1

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.p0.1.res \
      --random-seed $seed \
      -p 0 -t 1
if [ $? -gt 0 ]; then
    print_exit 1
fi

diff test/diff/$file.p0.1.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 2
fi

rm test/diff/$file.p0.1.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.p0.2.res \
      --random-seed $seed \
      -p 0 -t 2
if [ $? -gt 0 ]; then
    print_exit 3
fi

diff test/diff/$file.p0.2.res output_files/$file.p0.res
if [ $? -gt 0 ]; then
    print_exit 4
fi

rm test/diff/$file.p0.2.res

normal_exit
