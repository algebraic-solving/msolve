#!/usr/bin/env bash

# options --separating-element, --dense-route and --split on cyclic5-qq:
# the last variable is not separating, a separating linear form is found
# (0 4 3 2 1) and the multi-modular steps use the original system; the
# solutions split into 6 components, the largest one of degree 20 out of
# 70, so that the splitting is not used by default. Without the search,
# the changes of variable order fail and the linear form 1 2 3 4 5 is
# added: the parametrization differs. The other options do not change it.

file=cyclic5-qq

source test/diff/diff_source.sh

run_options() {
    local ex=$1 ref=$2
    shift 2
    $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.opt$ex.res \
          --random-seed $seed \
          -P 2 -d 4 -L 0 -l 2 "$@"
    if [ $? -gt 0 ]; then
        print_exit $ex
    fi

    diff test/diff/$file.opt$ex.res output_files/$ref
    if [ $? -gt 0 ]; then
        print_exit $((ex + 1))
    fi

    rm test/diff/$file.opt$ex.res
}

run_options 1 $file.P2.d4.sep0.res --separating-element 0 -t 1
run_options 3 $file.P2.d4.sep0.res --separating-element 0 -t 2
run_options 5 $file.P2.d4.res --separating-element 1 -t 1
run_options 7 $file.P2.d4.res --dense-route 1 -t 1
run_options 9 $file.P2.d4.res --dense-route 1 -t 2
run_options 11 $file.P2.d4.res --dense-route 2 -t 1
run_options 13 $file.P2.d4.res --split 1 -t 1
run_options 15 $file.P2.d4.res --split 2 -t 1
run_options 17 $file.P2.d4.res --split 2 -t 2

normal_exit
