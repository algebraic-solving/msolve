source test/diff/diff_source-param-fullnf.sh

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.81.res \
      --random-seed $seed \
      -P 2 -d 4 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 81
fi

diff test/diff/$file.81.res output_files/$file.P2.d4.res
if [ $? -gt 0 ]; then
    exit 82
fi

rm test/diff/$file.81.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.101.res \
      --random-seed $seed \
      -P 2 -d 4 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 101
fi

diff test/diff/$file.101.res output_files/$file.P2.d4.res
if [ $? -gt 0 ]; then
    exit 102
fi

rm test/diff/$file.101.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.121.res \
      --random-seed $seed \
      -P 2 -d 4 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 121
fi

diff test/diff/$file.121.res output_files/$file.P2.d4.res
if [ $? -gt 0 ]; then
    exit 122
fi

rm test/diff/$file.121.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.141.res \
      --random-seed $seed \
      -P 2 -d 4 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 141
fi

diff test/diff/$file.141.res output_files/$file.P2.d4.res
if [ $? -gt 0 ]; then
    exit 142
fi

rm test/diff/$file.141.res
