#!/usr/bin/env bash

# integer coefficients over GF(256) stand for their residues modulo 2,
# the expected result is the Groebner basis over GF(2)

file=random6-gf256-int

source test/diff/diff_source.sh

source test/diff/diff_source-gb-red.sh

normal_exit
