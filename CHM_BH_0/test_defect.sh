#!/usr/bin/bash

# 2022-08-02
#
# defect tester
#
# exemplary input:
# $1 = BH-12-12.txt (original BUTSON file)
# $2 = defect_BH_12_12.txt (calculated defects)
# $3 = 12 = q = BUTSON class
# $4 = N = number of tests
#
# $ ./test_defect.sh BH-12-12.txt defect_12_12.txt 12 100000 # run overnight...



MAX_LINE=$(wc -l $1 | awk NF=1);

for i in $(eval echo "{1..$4}"); do
    LINE_TO_READ=$((($RANDOM*$RANDOM)%$MAX_LINE))
    echo "checking defect in line: "$LINE_TO_READ
    ./check_defect.sh $1 $3 $LINE_TO_READ $2
done;

