#!/usr/bin/env bash

#test file for Issue #310: bivariate variant of
#univariate-cubic-close-roots.ms (Issue #358), pairing the same cubic
#in y (one exact rational root, two irrational roots extremely close
#together) with x=0, to check real root isolation determines the
#sign of every coordinate correctly in the genuinely multivariate
#code path too

file=issue-310-bivariate

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
