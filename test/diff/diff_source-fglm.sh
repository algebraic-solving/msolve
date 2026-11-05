$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.7.res \
      --random-seed $seed \
      -g 3 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    print_exit 7
fi

diff test/diff/$file.7.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 8
fi

rm test/diff/$file.7.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.27.res \
      --random-seed $seed \
      -g 3 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    print_exit 27
fi

diff test/diff/$file.27.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 28
fi

rm test/diff/$file.27.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.47.res \
      --random-seed $seed \
      -g 3 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    print_exit 47
fi

diff test/diff/$file.47.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 48
fi

rm test/diff/$file.47.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.67.res \
      --random-seed $seed \
      -g 3 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    print_exit 67
fi

diff test/diff/$file.67.res output_files/$file.g3.res
if [ $? -gt 0 ]; then
    print_exit 68
fi

rm test/diff/$file.67.res
