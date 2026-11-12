#!/usr/bin/env bash

# random dense quadratic system over GF(2), the expected result was
# computed with the linear algebra for primes < 2^8

file=random8-gf2

source test/diff/diff_source.sh

source test/diff/diff_source-gb-red.sh

normal_exit
