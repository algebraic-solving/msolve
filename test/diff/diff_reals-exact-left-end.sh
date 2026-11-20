#!/usr/bin/env bash

#test file for the real roots of x (2^256 x + 1) (3 2^256 x + 1): the
#isolating interval of the root -1/(3 2^256) starts at the root -2^-256
#of the same factor, and overlaps the root 0 of another factor, so it is
#halved until they are disjoint.

file=reals-exact-left-end

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
