$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.11.res \
      --random-seed $seed \
      -p 128 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 11
fi

diff test/diff/$file.11.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 12
fi

rm test/diff/$file.11.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.13.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 13
fi

diff test/diff/$file.13.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 14
fi

rm test/diff/$file.13.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.31.res \
      --random-seed $seed \
      -p 128 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 31
fi

diff test/diff/$file.31.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 32
fi

rm test/diff/$file.31.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.33.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 33
fi

diff test/diff/$file.33.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 34
fi

rm test/diff/$file.33.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.51.res \
      --random-seed $seed \
      -p 128 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 51
fi

diff test/diff/$file.51.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 52
fi

rm test/diff/$file.51.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.53.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 53
fi

diff test/diff/$file.53.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 54
fi

rm test/diff/$file.53..res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.71.res \
      --random-seed $seed \
      -p 128 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 71
fi

diff test/diff/$file.71.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 72
fi

rm test/diff/$file.71.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.73.res \
      --random-seed $seed \
      -p 1024 -L 0 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 73
fi

diff test/diff/$file.73.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 74
fi

rm test/diff/$file.73.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.91.res \
      --random-seed $seed \
      -p 128 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 91
fi

diff test/diff/$file.91.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 92
fi

rm test/diff/$file.91.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.93.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 93
fi

diff test/diff/$file.93.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 94
fi

rm test/diff/$file.93.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.111.res \
      --random-seed $seed \
      -p 128 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 111
fi

diff test/diff/$file.111.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 112
fi

rm test/diff/$file.111.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.113.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 113
fi

diff test/diff/$file.113.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 114
fi

rm test/diff/$file.113.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.131.res \
      --random-seed $seed \
      -p 128 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 131
fi

diff test/diff/$file.131.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 132
fi

rm test/diff/$file.131.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.133.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 133
fi

diff test/diff/$file.133.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 134
fi

rm test/diff/$file.133.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.151.res \
      --random-seed $seed \
      -p 128 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 151
fi

diff test/diff/$file.151.res output_files/$file.p128.res
if [ $? -gt 0 ]; then
    exit 152
fi

rm test/diff/$file.151.res

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.153.res \
      --random-seed $seed \
      -p 1024 -L 1 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 153
fi

diff test/diff/$file.153.res output_files/$file.p1024.res
if [ $? -gt 0 ]; then
    exit 154
fi

rm test/diff/$file.153.res
