#!/bin/bash

#test file for real root extraction with adaptative precision
#For this test file, one coordinate of one solutions is very close to 0
#msolve's real root extraction should be able to make the distinction

file=realroot1

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
