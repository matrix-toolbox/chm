function Y = SH_8_5_26(p)
% ------------------------------------------------------------------------------
% 2023-02-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% Symmetric 2-parameter family of CHM of order N = 8
% with generic d = 5 and #L = 26.
% It stems from BH(8, 2).
% Originally found as a special solution of "SH_8_3_142.m".
% -----------------------------------------------------------------------------
    a = exp(2j*pi* p(1));
    b = exp(2j*pi* p(2));

    Y = [
        1   1     1       1      1       1     1    1    ;
        1   a     b      -a     -b      -b/a  -1    b/a  ;
        1   b    -b*b/a  -b      b*b/a   b/a  -1   -b/a  ;
        1  -a    -b       a      b      -b/a  -1    b/a  ;
        1  -b     b*b/a   b     -b*b/a   b/a  -1   -b/a  ;
        1  -b/a   b/a    -b/a    b/a    -1     1   -1    ;
        1  -1    -1      -1     -1       1     1    1    ;
        1   b/a  -b/a     b/a   -b/a    -1     1   -1    ;
    ];

end

