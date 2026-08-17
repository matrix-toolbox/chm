function Y = SH_8_11_6A(p)
% ------------------------------------------------------------------------------
% 2023-02-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric 1-parameter family of CHM of order N = 8
% with generic d = 11 and #L = 6.
% It stems from BH(8, 2).
% It is a special solution of "SH_8_3_88A.m".
% ------------------------------------------------------------------------------

    a = exp(2j*pi * p(1));

    Y = [
        1   1   1   1   1   1   1   1  ;
        1   a  -1  -1   1  -1  -a   1  ;
        1  -1  -1   1  -1   1  -1   1  ;
        1  -1   1   1   1  -1  -1  -1  ;
        1   1  -1   1  -1  -1   1  -1  ;
        1  -1   1  -1  -1  -a'  1   a' ;
        1  -a  -1  -1   1   1   a  -1  ;
        1   1   1  -1  -1   a' -1  -a' ;
    ];

end

