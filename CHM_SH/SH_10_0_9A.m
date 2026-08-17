function Y = SH_10_0_9A
% -----------------------------------------------------------------------------
% 2023-01-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% Symmetric CHM of order 10 with d = 0 and #L = 9.
% It is a special solution of "SH_10_0_349A.m".
%------------------------------------------------------------------------------

     a = +11/16 - 3j*sqrt(15)/16;
     b = -1/4   - 1j*sqrt(15)/4;
     c = -7/8   - 1j*sqrt(15)/8;

    Y = [
        1  1  1  1  1  1  1  1  1  1  ;
        1  a  a  b  b' b  b' b' b' c  ;
        1  a  c  b' b' b' b' b  b  a  ;
        1  b  b' b  c' a' b  a' b  b' ;
        1  b' b' c' b  b  a' a' b  b  ;
        1  b  b' a' b  b  c' b  a' b' ;
        1  b' b' b  a' c' b  b  a' b  ;
        1  b' b  a' a' b  b  c' b  b' ;
        1  b' b  b  b  a' a' b  c' b' ;
        1  c  a  b' b  b' b  b' b' a  ;
    ];

end
