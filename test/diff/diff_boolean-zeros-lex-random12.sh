#!/usr/bin/env bash

# Groebner basis in the boolean ring w.r.t. the lexicographical order
# computed via the common zeros of the input polynomials (options
# --boolean-zeros --lex)

file=boolean-random12

source test/diff/diff_source.sh

for t in 1 2; do
    $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.lex.$t.res \
          --random-seed $seed \
          --boolean-zeros --lex -g 2 -t $t
    if [ $? -gt 0 ]; then
        print_exit $((2*t-1))
    fi

    diff test/diff/$file.lex.$t.res output_files/$file.lex.g2.res
    if [ $? -gt 0 ]; then
        print_exit $((2*t))
    fi

    rm test/diff/$file.lex.$t.res
done

normal_exit
