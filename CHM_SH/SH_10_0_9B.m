function Y = SH_10_0_9B
% -----------------------------------------------------------------------------
% 2023-01-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% Symmetric CHM of order 10 with d = 0 and #L = 9.
% It is a special solution of "SH_10_0_349B.m".
%------------------------------------------------------------------------------

    a = +17/32 + 1j*sqrt(735)/32;
    b = -1/4   + 1j*sqrt(15)/4;
    c = -7/8   - 1j*sqrt(15)/8;
    d = +11/16 - 1j*sqrt(135)/16;

    Y = [
        1  1    1    1      1      1    1    1   1       1    ;
        1  d    a    b      a      c    b^2  1   c       c    ;
        1  a    a    d      c      1    b^2  c   c       b    ;
        1  b    d    c      a*b/c  b    b    b   b'      b'   ;
        1  a    c    a*b/c  c^2    c    1    c   c/b     c    ;
        1  c    1    b      c      b'   1    c'  1       c'   ;
        1  b^2  b^2  b      1      1    b'   c'  c/a     1    ;
        1  1    c    b      c      c'   c'   b'  1       1    ;
        1  c    c    b'     c/b    1    c/a  1   b'*c/b  c/a  ;
        1  c    b    b'     c      c'   1    1   c/a     1    ;
    ];

end

