function Y = SH_9_4_30(p)
% ------------------------------------------------------------------------------
% 2023-02-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric 1-parameter family of CHM of order N = 9
% with generic d = 4 and #L = 30
% obtained from more general pattern: "SH_9_4_30_pattern.m".
% It stemms from BH(9, 6).
% ------------------------------------------------------------------------------

    a = exp(2j*pi * p(1));
    b = exp(2j*pi/3);
    d = b*a;
    c = b/a;
    Y = [
        1    1    1    1    1    1    1    1    1    ;
        1    b   -b'   a    c    -a   b'  -c    b'   ;
        1   -b'   b    b'   b'   b'  -b'   b'  -b'   ;
        1    a    b'   b    b'   b    d    1    c'   ;
        1    c    b'   b'   b    1    a'   b    d'   ;
        1   -a    b'   b    1    b   -d    b'  -c'   ;
        1    b'  -b'   d    a'  -d    b   -a'   b'   ;
        1   -c    b'   1    b    b'  -a'   b   -d'   ;
        1    b'  -b'   c'   d'  -c'   b'  -d'   b    ;
    ];

end

