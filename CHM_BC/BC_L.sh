#!/bin/bash
# 2021-05-14
# 2021-09-14
# 2022-04-18 full automation
# 2023-01-07 improved functionality
# 2026-09-13 output is named by the catalog convention: LH_N_d_L.data
#
# ************************************************************************
# *                                                                      *
# * THIS WORKS ONLY FOR N = 3 + 4*k FOR k = 0, 1, 2, ...                 *
# * One can try to find a solution for example for N = 13.               *
# * Users of Fractint should remember: "Higher orders take longer..." :) *
# *                                                                      *
# ************************************************************************
#
# Output is sorted in order to automatize the process of "unification".
# Given dimension N, prepare a circulant array template of order (N-1)/2;
# there should be total of (N+1)/2 unique equations, some of them must be
# removed manually...
#
#
# usage:
# $ ./BC_L.sh 11
# cat BC_L_DEBUG.txt
# octave >> H = LH_11;
#
# The generated M-file solves the system, then saves the phase matrix under
# the name the Catalog uses:  LH_<N>_<d>_<#Lambda>.data

N=$1
# WATCH OUT! For big N there might be a problem with leading zeros, however...
for j in $(seq -w 01 $(((N-1)/2))); do chars+=("c$j"); done;
#echo ${chars[@]} # DEBUG
#chars=(c01 c02 c03 c04 c05)

N="${#chars[@]}"
N2=$((1+2*N))
truncate -s 0 xxx1.txt
for j in $(seq 1 $N2)
do
        printf "1 " >> xxx1.txt
done
printf "\n" >> xxx1.txt

for j in $(seq 1 $N)
do
#        echo "${chars[@]}"
        printf "1 " >> xxx1.txt
        for ch in "${chars[@]}"
        do
                printf "%s %s" "$ch ${ch^^}" >> xxx1.txt
        done
        printf "\n" >> xxx1.txt
        printf "1 " >> xxx1.txt
        for ch in "${chars[@]}"
        do
                printf "%s %s" "${ch^^} $ch" >> xxx1.txt
        done
        printf "\n" >> xxx1.txt
     
        chars=("${chars[-1]}" "${chars[@]:0:${#chars[@]}-1}") # next shift
done

cp xxx1.txt BC_L_DEBUG.txt

T_FILE="xxx0.txt"
M_FILE="LH_"$N2".m"
truncate -s 0 $T_FILE

next_equ=0
readarray LSqr < xxx1.txt
k=${#LSqr[@]}
k=$((k-1))
#echo ""
#echo "array (core) size is "$k
#echo "there will be "$(((k+1)*k/2))" equations"
#echo ""


for r1 in $( eval echo {0..$((k-1))} )
do
	for r2 in $( eval echo {$((r1+1))..$((k))} )
	do
		next_equ=$((next_equ+1))
#		printf "z("$next_equ")="
		#echo ${LSqr[$r1]}" "${LSqr[$r2]}
		declare row1=(${LSqr[$r1]})
		declare row2=(${LSqr[$r2]})
		declare -a arr
		for j in $( eval echo {0..$k} )
		do
			declare e1=${row1[$j]};
			declare e2=${row2[$j]};
			if [[ $e1 =~ [1] ]] && [[ $e2 =~ [1] ]]; then
				arr[$j]=""
			fi
			if [[ $e1 =~ [1] ]] && [[ $e2 =~ [a-z] ]]; then
				arr[$j]=" +1/"$e2
			fi
			if [[ $e1 =~ [1] ]] && [[ $e2 =~ [A-Z] ]]; then
				arr[$j]=" +"$e2
			fi
			if [[ $e1 =~ [a-z] ]] && [[ $e2 =~ [1] ]]; then
				arr[$j]=" +"$e1
			fi
			if [[ $e1 =~ [A-Z] ]] && [[ $e2 =~ [1] ]]; then
				arr[$j]=" +1/"$e1
			fi
			if [[ $e1 =~ [a-z] ]] && [[ $e2 =~ [A-Z] ]]; then
				arr[$j]=" +"$e1"*"$e2
			fi
			if [[ $e1 =~ [a-z] ]] && [[ $e2 =~ [a-z] ]]; then
				arr[$j]=" +"$e1"/"$e2
			fi
			if [[ $e1 =~ [A-Z] ]] && [[ $e2 =~ [a-z] ]]; then
				arr[$j]=" +1"/$e1"/"$e2
			fi
			if [[ $e1 =~ [A-Z] ]] && [[ $e2 =~ [A-Z] ]]; then
				arr[$j]=" +"$e2"/"$e1
			fi
		done
	arr=( "${arr[@],,}" ) # to lowercase
	IFS=$'\n' sorted=($(sort <<<"${arr[*]}")); unset IFS
	printf "1" >> $T_FILE
	printf "%s" "${sorted[@]}" >> $T_FILE
#	printf "${e1}x${e2} "
        printf ";\n" >> $T_FILE
        done
done
	
#echo ""


sort -u $T_FILE > $M_FILE # remove identical lines

# add "y(j)=" at the begining of each line
truncate -s 0 $T_FILE
j=0
cat $M_FILE | while read -r line; do
        printf "\t" >> $T_FILE
        printf "%s" "y("$((j=j+1))")="$line >> $T_FILE
        printf "\n" >> $T_FILE
done


# replace characters with arguments: x(k)
for j in $(eval echo "{01..$N}")
do
	sed -i 's/c'$j'/x('$j')/g' $T_FILE
	# this produces some extra leading zeros; x(01) .. x(09), but they are transparent for Octave
	# we must keep leading zeros in: chars=(c01 c02 c03 ...) in the input
	# to avoid problems with order of matching patterns, eg. c1 and c11
done
sed -i 's/ //g' $T_FILE # remove spaces too

for j in $(eval echo "{01..$N}")
do
	sed -i 's/c'$j'/x('$j')/g' xxx1.txt
done

for j in $(eval echo "{01..$N}")
do
	sed -i 's/C'$j'/conj(x('$j'))/g' xxx1.txt
done

echo "" >> xxx1.txt # add one line for proper action of the next loop
truncate -s 0 xxx2.txt
cat xxx1.txt | while read -r line; do
        printf "\t\t" >> xxx2.txt
        printf "%s " $line";" >> xxx2.txt
        printf "\n" >> xxx2.txt
done



# prepare final form of M-file, eg. >> H = BC_L11;
echo "function L=LH_"$N2"()
% 2023-01-07
%
    clc
    N = "$N2";
    printf('Wait, calculating solution...');
    do
        [x info]=fsolve(@uc, exp(2j*pi*rand(1, (N-1)/2)));
        for j=1:64
            [x info]=fsolve(@uc, x);
        end
        s = norm(abs(x) - 1, 'fro');
        % unit modulus alone is not a Hadamard test -- the all-ones vector
        % passes it exactly and builds the all-ones matrix -- so the residual
        % of the orthogonality constraints has to be small as well
        r = norm(uc(x), 'fro');
        printf('.');
    until s < 1e-13 && r < 1e-13
    printf(' Solved!\n');
    L=[" > $M_FILE
cat xxx2.txt >> $M_FILE
echo "        ];
    A = mod(angle(L) + 2 * pi, 2 * pi) / 2 / pi;
    save(bc_name('LH', L, A), 'A');
end


function y=uc(x)" >> $M_FILE
cat $T_FILE >> $M_FILE
printf "%s\n\n" "end" >> $M_FILE

rm -rf xxx1.txt
rm -rf xxx2.txt
rm -rf $T_FILE
