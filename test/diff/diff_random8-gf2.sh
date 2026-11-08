#!/usr/bin/env bash

# random dense quadratic system over GF(2), the expected result was
# computed with the linear algebra for primes < 2^8

file=random8-gf2

source test/diff/diff_source.sh

for t in 1 4; do
    $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.$t.res \
          --random-seed $seed \
          -g 2 -t $t
    if [ $? -gt 0 ]; then
        print_exit $((2*t-1))
    fi

    diff test/diff/$file.$t.res output_files/$file.g2.res
    if [ $? -gt 0 ]; then
        print_exit $((2*t))
    fi

    rm test/diff/$file.$t.res
done

normal_exit
