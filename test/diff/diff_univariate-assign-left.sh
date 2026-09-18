#!/usr/bin/env bash

#test file for the assign_left bug fix (multivariate system whose
#rational univariate representation stresses usolve's real root
#isolation via assign_left)

file=univariate-assign-left

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
