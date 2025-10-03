$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -g 2 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 3
fi

diff test/diff/$file.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    exit 4
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -g 2 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 23
fi

diff test/diff/$file.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    exit 24
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -g 2 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 43
fi

diff test/diff/$file.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    exit 44
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -g 2 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 63
fi

diff test/diff/$file.res output_files/$file.g2.res
if [ $? -gt 0 ]; then
    exit 64
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 1 -g 2 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 31
fi

diff test/diff/$file.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 41
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 1 -g 2 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 231
fi

diff test/diff/$file.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 241
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 1 -g 2 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 431
fi

diff test/diff/$file.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 441
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 1 -g 2 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 631
fi

diff test/diff/$file.res output_files/$file.g2.e1.res
if [ $? -gt 0 ]; then
    exit 641
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 2 -g 2 -l 2 -t 1
if [ $? -gt 0 ]; then
    exit 32
fi

diff test/diff/$file.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 42
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 2 -g 2 -l 2 -t 2
if [ $? -gt 0 ]; then
    exit 232
fi

diff test/diff/$file.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 242
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 2 -g 2 -l 44 -t 1
if [ $? -gt 0 ]; then
    exit 432
fi

diff test/diff/$file.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 442
fi

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res \
      -e 2 -g 2 -l 44 -t 2
if [ $? -gt 0 ]; then
    exit 632
fi

diff test/diff/$file.res output_files/$file.g2.e2.res
if [ $? -gt 0 ]; then
    exit 642
fi

rm test/diff/$file.res
