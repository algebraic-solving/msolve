seed=${SEED:-$EPOCHSECONDS}

if [ -t 2 ]; then
    colseed='\033[0;96m' # High Intensity Light blue.
    colexit='\033[0;93m' # High Intensity Yellow.
    std='\033[0m'
else
    colseed=
    colexit=
    std=
fi

echo -e "${colseed}SEED${std}: $seed"

# each diff_example.sh is built by running msolve on $file.ms
# with options -L 0 -l 2 -t 1
# if the execution fails, exit 1
# then compare the output with the expected one, if different, exit 2
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
