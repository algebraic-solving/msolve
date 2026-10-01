#!/usr/bin/env bash

# Probabilistic sparse/dense linear algebra (-l 42) used to lose rows on
# small matrices since the random multipliers were not uniformly chosen
# in [1, p-1]. Then msolve wrote an empty output file while exiting with
# status 0. The failure depends on the random seed, thus we run several
# consecutive seeds. Line 5 of the output (the random linear form) depends
# on the seed and is not compared.

file=nonradical-radicalshape-31

source test/diff/diff_source.sh

sed 5d output_files/$file.res > test/diff/$file.prob-la.ref

for i in $(seq 0 39); do
    for t in 1 2; do
        $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.prob-la.res \
            --random-seed $((seed + i)) -P 1 -l 42 -t $t
        if [ $? -gt 0 ]; then
            print_exit 1
        fi

        sed 5d test/diff/$file.prob-la.res | \
            diff - test/diff/$file.prob-la.ref
        if [ $? -gt 0 ]; then
            print_exit 2
        fi

        rm test/diff/$file.prob-la.res
    done
done

rm test/diff/$file.prob-la.ref

normal_exit
