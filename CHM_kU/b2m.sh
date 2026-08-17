#!/usr/bin/bash

# 2022-02-04
# 2022-07-28
#
# read n^th line from BH-NN-qq.txt and convert it into M-FILE
#
# usage:
# $ ./b2m BH-NN-qq.txt q line_number # where q = BUTSON parameter
#
# -----------------------------------------------------------------------------
# possible issues:
# * before first run check the path to bash: #!/usr/home/bash
# * chmod +x *.sh
# -----------------------------------------------------------------------------

M_FILE=TBH.m # <--- temporary file with BH

N=$(echo $1 | awk -F[--] '{print $2}') # get dimension from file name (always between "-"; eg. -12-
dots=$(for i in $(eval echo "{1..$N}"); do echo -n "."; done)
nth_line=$3
echo 'function B=TBH' > $M_FILE
##echo "warning('off', 'all');" >> $M_FILE
echo 'B=exp(2j*pi*[' >> $M_FILE

# for fucks sake... semicolon is temporarily replaced with #, which is replaced back to semicolon at the end, because semicolon is used to encode 11 (eleven)
sed "${nth_line}q;d" $1 | cut -c 7- | sed 's/\('$dots'\)/\1#\n/g' | sed 's/\(.\)/\1 /g' | sed 's/:/10/g' | sed 's/;/11/g' | sed 's/</12/g' | sed 's/=/13/g' | sed 's/>/14/g' | sed 's/?/15/g' | sed 's/@/16/g' | sed 's/#/;/g' >> $M_FILE

echo "]/$2)/sqrt($N);" >> $M_FILE # all Butsons are normalized to UNITARY!
#echo "disp(sprintf('%g\t%g', SL(R(B)), SL(T(B))))" >> $M_FILE
#echo 'defect=[ud(B,"S", 1e-8) ud(B, "T") ud(B, "R")]' >> $M_FILE
echo 'end' >> $M_FILE
