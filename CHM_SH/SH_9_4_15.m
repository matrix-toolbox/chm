function Y = SH_9_4_15(p)
% ------------------------------------------------------------------------------
% 2023-02-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric 1-parameter family of CHM of order N = 9
% with generic d = 9 and #L = 15.
% It stems from BH(9, 3).
% It is a special solution of "SH_9_0_89.m".
% ------------------------------------------------------------------------------

    w = exp(2j*pi/3);
    a = exp(2j*pi * p(1));
    b = a * w^2;

    Y = [
        1  1   1        1       1       1         1   1          1        ;
        1  1   w        w       w'      w         1   w'         w'       ;
        1  w   b^2      w*b^2   b       b^2/w     w'  a^2/b      b^3/a^2  ;
        1  w   w*b^2    a^2     a^2/b   a^2/w'    w'  a          b        ;
        1  w'  b        a^2/b   a^2     a         w   a^2/w'     w*b^2    ;
        1  w   b^2/w    a^2/w'  a       a^4/b^2   w'  a^4/b^2/a  a^2/b    ;
        1  1   w'       w'      w       w'        1   w          w        ;
        1  w'  a^2/b    a       a^2/w'  a^4/b^2/a w   a^4/b^2    b^2/w    ;
        1  w'  b^3/a^2  b       w*b^2   a^2/b     w   b^2/w      b^2      ;
    ];

end

