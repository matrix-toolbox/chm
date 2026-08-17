#!/usr/bin/bash

# 2022-07-29
#
# count number of zeros in defect_BH_NN_qq.txt file
# = number of isolated BH matrices
#
# usage:
# $ ./count_isolated.sh BH-NN-qq.txt


grep -w -c 0 $1
