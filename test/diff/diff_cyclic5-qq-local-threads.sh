#!/usr/bin/env bash

# Groebner basis over the rationals (-g 2) with several threads applying
# the tracer to each prime (--local-threads), which must divide the number
# of threads: otherwise msolve warns and uses one thread per prime.

file=cyclic5-qq

source test/diff/diff_source.sh

# 2 primes in parallel, 2 threads each
$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.lt1.res \
      --random-seed $seed \
      -g 2 -t 4 --local-threads 2
if [ $? -gt 0 ]; then
    print_exit 1
fi

diff test/diff/$file.lt1.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    print_exit 2
fi

rm test/diff/$file.lt1.res

# 1 prime at a time, 4 threads
$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.lt2.res \
      --random-seed $seed \
      -g 2 -t 4 --local-threads 4
if [ $? -gt 0 ]; then
    print_exit 3
fi

diff test/diff/$file.lt2.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    print_exit 4
fi

rm test/diff/$file.lt2.res

# 3 does not divide 4: warning, one thread per prime
$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.lt3.res \
      --random-seed $seed \
      -g 2 -t 4 --local-threads 3 2> test/diff/$file.lt3.err
if [ $? -gt 0 ]; then
    print_exit 5
fi

grep -q "do not divide" test/diff/$file.lt3.err
if [ $? -gt 0 ]; then
    print_exit 6
fi

diff test/diff/$file.lt3.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    print_exit 7
fi

rm test/diff/$file.lt3.res test/diff/$file.lt3.err

normal_exit
