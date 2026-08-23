#!/usr/bin/env bash

file=reject-malformed

source test/diff/diff_source.sh

expect_fail() {
    local in=$1
    local excode=$2
    $(pwd)/msolve -g 2 -t 1 --random-seed $seed \
        -f input_files/$in -o test/diff/$file.out
    if [ $? -eq 0 ]; then
        print_exit $excode
    fi
}

# Parentheses, post-monomial division, a leading "*", and a coefficient
# that does not fit in a machine word must all be rejected.
expect_fail bad-paren.ms 1
expect_fail bad-div.ms 3
expect_fail bad-star.ms 5
expect_fail bad-overflow-ff.ms 7

# Documented leading-rational form must still parse.
$(pwd)/msolve -g 2 -t 1 --random-seed $seed \
    -f input_files/good-leading-rat.ms -o test/diff/$file.good.res
if [ $? -gt 0 ]; then
    print_exit 9
fi
diff_gb_output test/diff/$file.good.res output_files/good-leading-rat.g2.res
if [ $? -gt 0 ]; then
    print_exit 10
fi
rm -f test/diff/$file.out test/diff/$file.good.res

normal_exit
