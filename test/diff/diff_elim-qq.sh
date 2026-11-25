#!/usr/bin/env bash

file=elim-qq

source test/diff/diff_source.sh

source test/diff/diff_source-gb-elim-red.sh

# tiny initial hash table (-s 1): it must be enlarged while learning the trace
$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.1.res \
      --random-seed $seed \
      -e 1 -g 2 -s 1 -t 1
if [ $? -gt 0 ]; then
    print_exit 1
fi

diff test/diff/$file.1.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    print_exit 2
fi

rm test/diff/$file.1.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.21.res \
      --random-seed $seed \
      -e 1 -g 2 -s 1 -t 2
if [ $? -gt 0 ]; then
    print_exit 21
fi

diff test/diff/$file.21.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    print_exit 22
fi

rm test/diff/$file.21.res

normal_exit
