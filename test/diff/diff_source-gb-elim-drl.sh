# each diff_elim-drl-*.sh runs msolve on $file.ms with an elimination
# order, the input has more equations than variables, thus by default
# a DRL Groebner basis is computed first (option --elim-drl 2) unless
# only one variable is eliminated
# if the execution fails, print_exit with an odd code
# then compare the output with the expected one, if different, print_exit
# with the following even code
run_elim_drl() {
    local code=$1
    local elim=$2
    shift 2
    $(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.$code.res \
        --random-seed $seed \
        -e $elim -g 2 "$@"
    if [ $? -gt 0 ]; then
        print_exit $code
    fi

    diff test/diff/$file.$code.res output_files/$file.g2.e$elim.res
    if [ $? -gt 0 ]; then
        print_exit $((code+1))
    fi

    rm test/diff/$file.$code.res
}

run_elim_drl 1 4 -t 1
run_elim_drl 3 4 -t 2
run_elim_drl 5 4 --elim-drl 1 -t 1
run_elim_drl 7 5 -t 1
run_elim_drl 9 5 -t 2
run_elim_drl 11 5 --elim-drl 1 -l 44 -t 2
run_elim_drl 13 1 -t 1
run_elim_drl 15 1 --elim-drl 1 -t 2
