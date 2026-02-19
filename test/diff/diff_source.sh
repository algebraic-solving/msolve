seed=${SEED:-$(date +%s)}

setup_output() {

    # already initialized? -> do nothing
    [ -n "${_msolve_test_output:-}" ] && return
    _msolve_test_output=1

    # defaults (safe for logs / CI / pipes)
    isatty=0
    use_color=0
    use_cr=0

    # No colored output (default)
    colseed=""
    colexit=""
    std=""

    # Last character after normal exit
    lastcharexit="\n"

    # Detect terminal on stderr (diagnostics channel)
    if [ -t 2 ] && [ -n "$TERM" ] && [ "$TERM" != dumb ]; then
        isatty=1
        use_color=1
        use_cr=1
    fi

    # Colored output
    if [ "$use_color" -eq 1 ]; then
	# colseed='\033[0;96m' # High Intensity Light blue.
	colseed=$(tput setaf 14) # High Intensity Light blue.
	# colexit='\033[0;93m' # High Intensity Yellow.
	colexit=$(tput setaf 11) # High Intensity Yellow.
	# std='\033[0m'
	std=$(tput sgr0)
    fi

    # Carriage return afer normal exit
    if [ "$use_cr" -eq 1 ]; then
	lastcharexit="\r"
    fi
}

# print the seed in color
print_seed() {
    setup_output
    printf "${beforeseed}${colseed}SEED${std}: %-10d " $seed >&2
}

print_seed

# print a last character depending on the output stream
normal_exit() {
    printf "$lastcharexit" >&2
}


# print the exit code in color
print_exit() {
    local excode=$1
    printf "${colexit}EXIT${std}: $excode\n" >&2
    exit "$excode"
}

# each diff_example.sh is built by running msolve on $file.ms
# with options -L 0 -l 2 -t 1
# if the execution fails, print_exit 1
# then compare the output with the expected one, if different, print_exit 2
# repeat for other execution parameters with exit 3 and 4 and so on, until at most 19 and 20

# repeat all these tests changing -L 0 into -L 1, -l 2 into -l 44 and -t 1 into -t 2
# increase all exit codes by
# 20 for -t 2
# 40 for -l 44
# thus 60 for -l 44 and -t 2
# 80 for -L 1
# thus 100 for -L 1 and -t 2
# thus 120 for -L 1 and -l 44
# thus 140 for -L 1, -l 44 and -t 2
