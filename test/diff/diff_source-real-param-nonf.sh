$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.1.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.1.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 2
fi

rm test/diff/$file.1.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.3.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.3.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 4
fi

rm test/diff/$file.3.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.21.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 21
fi

diff test/diff/$file.21.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 22
fi

rm test/diff/$file.21.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.23.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 23
fi

diff test/diff/$file.23.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 24
fi

rm test/diff/$file.23.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.41.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 41
fi

diff test/diff/$file.41.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 42
fi

rm test/diff/$file.41.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.43.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 43
fi

diff test/diff/$file.43.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 44
fi

rm test/diff/$file.43.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.61.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 61
fi

diff test/diff/$file.61.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 62
fi

rm test/diff/$file.61.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.63.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 63
fi

diff test/diff/$file.63.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 64
fi

rm test/diff/$file.63.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.81.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 81
fi

diff test/diff/$file.81.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 82
fi

rm test/diff/$file.81.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.83.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 83
fi

diff test/diff/$file.83.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 84
fi

rm test/diff/$file.83.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.101.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 101
fi

diff test/diff/$file.101.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 102
fi

rm test/diff/$file.101.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.103.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 103
fi

diff test/diff/$file.103.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 104
fi

rm test/diff/$file.103.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.121.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 121
fi

diff test/diff/$file.121.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 122
fi

rm test/diff/$file.121.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.123.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 123
fi

diff test/diff/$file.123.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 124
fi

rm test/diff/$file.123.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.141.res \
      --random-seed $seed \
      -P 1 -d 0 -p 128 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 141
fi

diff test/diff/$file.141.res output_files/$file.P1.d0.p128.res
if [ $? -gt 0 ]; then
    exit 142
fi

rm test/diff/$file.141.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.143.res \
      --random-seed $seed \
      -P 1 -d 0 -p 1024 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 143
fi

diff test/diff/$file.143.res output_files/$file.P1.d0.p1024.res
if [ $? -gt 0 ]; then
    exit 144
fi

rm test/diff/$file.143.res
