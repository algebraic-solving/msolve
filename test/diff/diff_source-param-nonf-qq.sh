source test/diff/diff_source-param-nonf.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 2 -d 0 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 81
fi

diff test/diff/$file.res output_files/$file.P2.d0.res
if [ $? -gt 0 ]; then
    exit 82
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 2 -d 0 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 101
fi

diff test/diff/$file.res output_files/$file.P2.d0.res
if [ $? -gt 0 ]; then
    exit 102
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 2 -d 0 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 121
fi

diff test/diff/$file.res output_files/$file.P2.d0.res
if [ $? -gt 0 ]; then
    exit 122
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 2 -d 0 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 141
fi

diff test/diff/$file.res output_files/$file.P2.d0.res
if [ $? -gt 0 ]; then
    exit 142
fi

rm test/diff/$file.res
