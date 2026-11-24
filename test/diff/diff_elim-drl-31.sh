#!/usr/bin/env bash

file=elim-drl-31

source test/diff/diff_source.sh

source test/diff/diff_source-gb-elim-drl.sh

# without DRL Groebner basis first
run_elim_drl 17 4 --elim-drl 0 -t 1

normal_exit
