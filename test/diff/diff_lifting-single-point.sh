#!/usr/bin/env bash

#test file for a system with a single solution, (0, 0, -1/3, 0, 2): the
#eliminating polynomial has degree 1 and its p-adic lifting must reduce
#the last coordinate t modulo it.

file=lifting-single-point

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
