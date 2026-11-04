#!/usr/bin/env bash

file=lexgb-nonshape-31

# reduced Groebner basis for the lex order by FGLM (-g 3): ideal not in shape position (x^2-1, y^2-1, z^2-1)
# exit codes: 1 (execution), 2 (diff), shifted by 20 for -t 2, 40 for -l 44,
# (see test/diff/diff_source.sh)

source test/diff/diff_source.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.1.res \
      -g 3 -l 2 -t 1
if [ $? -gt 0 ]; then
    print_exit 1
fi

diff test/diff/$file.1.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 2
fi

rm test/diff/$file.1.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.21.res \
      -g 3 -l 2 -t 2
if [ $? -gt 0 ]; then
    print_exit 21
fi

diff test/diff/$file.21.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 22
fi

rm test/diff/$file.21.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.41.res \
      -g 3 -l 44 -t 1
if [ $? -gt 0 ]; then
    print_exit 41
fi

diff test/diff/$file.41.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 42
fi

rm test/diff/$file.41.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.61.res \
      -g 3 -l 44 -t 2
if [ $? -gt 0 ]; then
    print_exit 61
fi

diff test/diff/$file.61.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 62
fi

rm test/diff/$file.61.res

normal_exit
