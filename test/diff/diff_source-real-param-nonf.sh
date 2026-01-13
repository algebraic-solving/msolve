$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 2
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 4
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 21
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 22
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 23
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 24
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 41
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 42
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 43
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 44
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 61
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 62
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 63
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 64
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 81
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 82
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 83
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 84
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 101
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 102
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 103
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 104
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 121
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 122
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 123
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 124
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 141
fi

diff test/diff/$file.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 142
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 143
fi

diff test/diff/$file.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 144
fi

rm test/diff/$file.res
