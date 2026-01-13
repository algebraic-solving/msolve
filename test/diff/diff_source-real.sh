$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 11
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 12
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 13
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 14
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 31
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 32
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 33
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 34
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 51
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 52
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 53
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 54
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 71
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 72
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 73
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 74
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 91
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 92
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 93
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 94
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 111
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 112
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 113
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 114
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 131
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 132
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 133
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 134
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 128 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 151
fi

diff test/diff/$file.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 152
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 153
fi

diff test/diff/$file.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 154
fi

rm test/diff/$file.res
