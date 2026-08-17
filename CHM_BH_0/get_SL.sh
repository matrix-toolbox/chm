#!/usr/bin/bash

# 2022-07-22
# usage:
# $ ./get_SL.sh SPLIT_FILEax RESULT_xx.txt
#
# THIS METHOD IS SLOW AND SHOULD NOT BE USED! USE "BHrc.ELF" INSTEAD!

for j in {1..100000}; do ./b2m.sh $1 4 $j; octave TBH.m; done > $2

