#!/usr/bin/env bash

# Groebner basis in the Boolean ring over GF(2) for an input file over
# GF(2) via the option --boolean, eliminating the first three variables

file=boolean-elim10

source test/diff/diff_source.sh

for t in 1 2; do
    $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.$t.res \
          --random-seed $seed \
          --boolean -e 3 -g 2 -t $t
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
