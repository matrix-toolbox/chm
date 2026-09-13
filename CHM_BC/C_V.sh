#!/bin/bash
# 2021-05-17
# 2022-04-17
# 2023-01-07 improved functionality
# 2026-09-13 output is named by the catalog convention: VH_N_d_L.data
#
# ***************************************************
# *                                                 *
# * THIS WORKS FOR _ANY_ DIMENSION N = 6, 7, 8, ... *
# *                                                 *
# ***************************************************

# prepare ready-to-use M-file to calculate UC from undephased circulant pattern of size N
# usage:
# $ ./C_V.sh 6
# octave >> H = VH_6;
#
# The generated M-file solves the system, then saves the phase matrix under
# the name the Catalog uses:  VH_<N>_<d>_<#Lambda>.data


N=$1
M_FILE="VH_"$N".m"
truncate -s 0 $M_FILE


# dollar and single quotes allows inserting escape characters: uc=$'\n'
# then must be printed in double quotes: echo "$uc"
uc=$'function y=uc(x)\n'
for k in $(seq 1 $((N-1)))
do
	uc="$uc y("$k")="
	for j in $(seq 1 $N)
	do
		t1="x("$j")/x("$(((j+k) % N))")"
		t1=$(echo $t1 | sed "s/(0)/("$N")/g")
	        uc="$uc $t1"
		if [[ $j -le $N-1 ]]
		then
			t2="+"
			uc="$uc $t2"
		fi

	done
	t3=$';\n'
	uc="$uc $t3"
	unset t1
	unset t2
	unset t3
done
uc="$uc""end"


Hm=""
for k in $(seq 1 $((N)))
do
	for j in $(seq 1 $N)
	do
		h1="x("$(((j-k+1+N) % N))")"
		h1=$(echo $h1 | sed "s/(0)/("$N")/g")
		Hm="$Hm $h1"
	done
	t3=$';\n'
	Hm="$Hm $t3"
	unset t3	
done


# prepare final form of M-file
echo "function V=VH_"$N"()
% 2023-01-07
    clc
    N = "$N";
    printf('Wait, calculating solution...');
    do
        [x info]=fsolve(@uc, exp(2j*pi*rand(1, N)));
        for j=1:64
            [x info]=fsolve(@uc, x);
        end
        x = x / abs(x(1));
        s = norm(abs(x) - 1, 'fro');
        printf('.');
    until s < 1e-13
    printf(' Solved!\n');
    K=[" > $M_FILE
echo "$Hm" >> $M_FILE
echo "        ];
    V = (diag(1 ./ K(:, 1)) * K);
    V = V * diag(1 ./ V(1, :));
    A = mod(angle(V) + 2 * pi, 2 * pi) / 2 / pi;
    save(bc_name('VH', V, A), 'A');
end" >> $M_FILE
echo $'\n'"$uc"$'\n\n' >> $M_FILE


