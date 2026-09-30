#!/usr/bin/env bash

# With this seed and glibc's rand(), the Hankel matrix used to compute the
# parametrizations has a singular sub-block of size dim - 1, so that the
# computation must fall back to another method.

file=cyclic5-16

SEED=1790795968

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/hankel-degenerate-16.res \
      --random-seed $seed \
      -P 2 -d 4 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    print_exit 1
fi

diff test/diff/hankel-degenerate-16.res output_files/$file.P2.d4.res
if [ $? -gt 0 ]; then
    print_exit 2
fi

rm test/diff/hankel-degenerate-16.res

normal_exit
