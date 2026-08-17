#!/usr/bin/bash
# 2022-08-03

# usage:
# $ ./get_lines_with_zero_defect.sh defect_BH_12_8.txt

grep -n --color=auto "^0" $1


# manual search for zeros:
#   n++: \r\n0\r\n
#   vim: /\n0
#
# ALWAYS keep NEWLINE at the end of the file!
#
#
#


