#!/usr/bin/env bash

#test file for real root isolation sign determination: one coordinate
#(y4) is a root extremely close to 0, msolve's real root isolation
#should still determine its sign correctly

file=univariate-sgn

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
