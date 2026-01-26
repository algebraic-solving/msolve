$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.5.res \
      --random-seed $seed \
      -e 1 -g 2 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 5
fi

diff test/diff/$file.5.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 6
fi

rm test/diff/$file.5.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.7.res \
      --random-seed $seed \
      -e 2 -g 2 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 7
fi

diff test/diff/$file.7.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 8
fi

rm test/diff/$file.7.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.25.res \
      --random-seed $seed \
      -e 1 -g 2 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 25
fi

diff test/diff/$file.25.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 26
fi

rm test/diff/$file.25.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.27.res \
      --random-seed $seed \
      -e 2 -g 2 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 27
fi

diff test/diff/$file.27.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 28
fi

rm test/diff/$file.27.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.45.res \
      --random-seed $seed \
      -e 1 -g 2 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 45
fi

diff test/diff/$file.45.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 46
fi

rm test/diff/$file.45.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.47.res \
      --random-seed $seed \
      -e 2 -g 2 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 47
fi

diff test/diff/$file.47.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 48
fi

rm test/diff/$file.47.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.65.res \
      --random-seed $seed \
      -e 1 -g 2 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 65
fi

diff test/diff/$file.65.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 66
fi

rm test/diff/$file.65.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.67.res \
      --random-seed $seed \
      -e 2 -g 2 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 67
fi

diff test/diff/$file.67.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 68
fi

rm test/diff/$file.67.res


