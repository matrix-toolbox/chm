function Y = SH_10_9_20(p)
% ------------------------------------------------------------------------------
% 2023-02-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric 1-parameter family of CHM of order N = 10
% with generic d = 9 and #L = 20.
% It stemms from BH(10, 4).
% It is a special solution of "SH_10_0_349B.m".
% ------------------------------------------------------------------------------

    a = exp(2j*pi* p(1));

    Y = [
        1    1    1    1    1    1    1    1    1    1    ;
        1   -1    a    I   -I    a   -1   -a    1   -a    ;
        1    a    1    a'   I   -1   -a   -1   -a'  -I    ;
        1    I    a'  -1   -I   -a'   1    a'  -1   -a'   ;
        1   -I    I   -I   -1    I   -I    I   -I    I    ;
        1    a   -1   -a'   I   -I   -a    1    a'  -1    ;
        1   -1   -a    1   -I   -a   -1    a    I    a    ;
        1   -a   -1    a'   I    1    a   -I   -a'  -1    ;
        1    1   -a'  -1   -I    a'   I   -a'  -1    a'   ;
        1   -a   -I   -a'   I   -1    a   -1    a'   1    ;
    ];
end

