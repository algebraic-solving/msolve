#!/bin/bash

file=bug_68

diff test/diff/$file.res output_files/$file.res -l 2
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res -l 44
if [ $? -gt 0 ]; then
    exit 1
fi

rm test/diff/$file.res
