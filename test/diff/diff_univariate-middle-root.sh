#!/usr/bin/env bash

#test file for real root isolation: the quintic has the exact root
#171/256, at the middle of an interval of the bisection, and two roots
#closer than 2^-18 to each other. At the interval whose middle is
#171/256, the right half has no root and is skipped by the bisection.

file=univariate-middle-root

source test/diff/diff_source.sh

source test/diff/diff_source-real.sh

normal_exit
