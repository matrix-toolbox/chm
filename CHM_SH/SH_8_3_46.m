function Y = SH_8_3_46(p)
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric 2-parameter family of CHM of order N = 8
% with generic d = 3 and #L = 46.
% It stems from BH(8, 2).
% It was originally found as a solution of "SH_8_3_142.m".
% ------------------------------------------------------------------------------

    a = exp(2j*pi* p(1));
    b = exp(2j*pi* p(2));

    c = -(1+a^2+2*b)/(2+b+b/a^2);

    Y = [
        1   1     1     1     1         1      1         1   ;
        1   c/a  -a    -1     c         a     -c/a      -c   ;
        1  -a    -a^2   a     a^2      -a     -1         a   ;
        1  -1     a    -b/a   b        -b      b/a      -a   ;
        1   c     a^2   b     c*b       b      c*b/a^2   c   ;
        1   a    -a    -b     b         b/a   -b/a      -1   ;
        1  -c/a  -1     b/a   c*b/a^2  -b/a   -c*b/a^2   c/a ;
        1  -c     a    -a     c        -1      c/a      -c/a ;
    ];

end

