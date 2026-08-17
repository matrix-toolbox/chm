function Y = HH_12_11_6
% ------------------------------------------------------------------------------
% 2023-01-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian BH(12, 6) with d = 11 and #L = 6.
% Found as a solution of "HH_12_5_1902.m" and defined as a separate pattern.
% ------------------------------------------------------------------------------

    b = exp(8j*pi/6);
    Y = [
        1   1   1   1   1   1   1   1   1   1   1   1  ;
        1   1   b   b   b'  b'  b   b'  1   1   b   b' ;
        1   b'  1   b   b'  b   b'  1   b   b   1   b' ;
        1   b'  b'  1   1   b'  b'  b   1   b   b   b  ;
        1   b   b   1   1   b   1   b'  b'  b   b'  b' ;
        1   b   b'  b   b'  1   1   1   b'  b'  b   b  ;
        1   b'  b   b   1   1  -1  -1  -1  -b  -b  -b' ;
        1   b   1   b'  b   1  -1  -1  -b  -1  -b' -b  ;
        1   1   b'  1   b   b  -1  -b' -1  -b  -1  -b  ;
        1   1   b'  b'  b'  b  -b' -1  -b' -1  -b  -b' ;
        1   b'  1   b'  b   b' -b' -b  -1  -b' -1  -b' ;
        1   b   b   b'  b   b' -b  -b' -b' -b  -b  -1  ;
    ];
    
end
