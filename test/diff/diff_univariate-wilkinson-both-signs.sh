#!/usr/bin/env bash

#test file for real root isolation: x (prod_{0 < |i| <= 30} (x - i) - 1)
#has 61 real roots, of both signs and of absolute value up to 30. Parts
#of the bisection are left to the threads, and the roots they find are
#put back in order.

file=univariate-wilkinson-both-signs

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
