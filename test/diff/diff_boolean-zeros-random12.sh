#!/usr/bin/env bash

# Groebner basis in the boolean ring computed via the common zeros of the
# input polynomials (option --boolean-zeros), it has to be the same as
# the one computed by F4

file=boolean-random12

source test/diff/diff_source.sh

for t in 1 2; do
    $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.zeros.$t.res \
          --random-seed $seed \
          --boolean-zeros -g 2 -t $t
    if [ $? -gt 0 ]; then
        print_exit $((2*t-1))
    fi

    diff test/diff/$file.zeros.$t.res output_files/$file.g2.res
    if [ $? -gt 0 ]; then
        print_exit $((2*t))
    fi

    rm test/diff/$file.zeros.$t.res
done

normal_exit
