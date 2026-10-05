#!/usr/bin/env bash

file=elim-full-qq

source test/diff/diff_source.sh

# full Groebner basis w.r.t. an elimination order over the rationals,
# not only the basis of the elimination ideal
for e in 1 2; do
    for t in 1 2; do
        for l in 2 44; do
            $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.$e.$t.$l.res \
                  --random-seed $seed \
                  -e $e --elim-full-basis -g 2 -l $l -t $t
            if [ $? -gt 0 ]; then
                print_exit 1
            fi

            diff test/diff/$file.$e.$t.$l.res output_files/$file.g2.e$e.res
            if [ $? -gt 0 ]; then
                print_exit 2
            fi

            rm test/diff/$file.$e.$t.$l.res
        done
    done
done

normal_exit
