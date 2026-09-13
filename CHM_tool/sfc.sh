#!/usr/bin/octave -qf

# 2022-08-20
#
# Quick SFC-CHECK of data saved by "sinkhorn*.m" procedure.
#
# usage:
# $ chmod +x sfc.sh
#
# $ for j in $(ls Y10_4*T*.dat); do ./sfc.sh $j; done >> FILE_NAME.dat
# $ tail -f FILE_NAME.dat
#
# or:
# $ ./sfc.sh Y7_0_L43_20220902T230104.dat

arg_list = argv();

addpath ../CHM

fn = arg_list{1}
load(fn);
sfc(exp(2j*pi*A), "VERBOSE")
