#!/usr/bin/env bash

# integer coefficients over GF(256) stand for their residues modulo 2,
# the expected result is the Groebner basis over GF(2)

file=random6-gf256-int

source test/diff/diff_source.sh

for t in 1 2; do
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
