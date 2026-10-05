#!/usr/bin/env bash

source test/diff/diff_source.sh

# the output file should use the same line endings as the input file:
# for each pair of inputs differing only by their line endings, the output
# for the dos input must be the output for the unix input with "\r\n"

prefix=input_files/line_endings

check_line_endings() {
    local dos=$1
    local unix=$2
    local excode=$3

    $(pwd)/msolve -f $prefix/$dos.ms -o test/diff/$dos.res \
          --random-seed $seed $4
    if [ $? -gt 0 ]; then
        print_exit $excode
    fi

    $(pwd)/msolve -f $prefix/$unix.ms -o test/diff/$unix.res \
          --random-seed $seed $4
    if [ $? -gt 0 ]; then
        print_exit $(($excode+1))
    fi

    if grep -q $'\r' test/diff/$unix.res; then
        print_exit $(($excode+2))
    fi

    awk '{ printf "%s\r\n", $0 }' test/diff/$unix.res | \
        cmp -s - test/diff/$dos.res
    if [ $? -gt 0 ]; then
        print_exit $(($excode+3))
    fi

    rm test/diff/$dos.res test/diff/$unix.res
}

check_line_endings in1_dos in1_unix 1
check_line_endings in2_dos_noeol in2_unix 5
check_line_endings in3_dos in3_unix 9
check_line_endings in4_dos in4_unix 13
check_line_endings in1_dos in1_unix 17 "-g 2"
check_line_endings in1_dos in1_unix 21 "-P 1"

normal_exit
