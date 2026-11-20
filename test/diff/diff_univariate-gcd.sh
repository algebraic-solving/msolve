#!/usr/bin/env bash

#test file for univariate input: two polynomials with rational
#coefficients whose gcd (x^2-2)^2*(x-1) is not squarefree. The solutions
#are given by its squarefree part, without Groebner bases.

file=univariate-gcd

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
