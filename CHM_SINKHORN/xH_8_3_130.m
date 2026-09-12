function Y = xH_8_3_130
% ------------------------------------------------------------------------------
% 2022-01-19
% 2025-06-08
%
% complex Hadamard matrix of order 8
%
% c1 and c2 are fixed exotic roots of unity
% a, c, d, e are numeric constants with no closed form identified yet
% -- not symmetrizable
% ------------------------------------------------------------------------------

    c1 = exp(23j*pi/41);
    c2 = exp(17j*pi/53);

    a =  -0.70388420081681013 + 1j*0.71031474139319428;
    c =   0.80929591242166365 - 1j*0.58740116286707067;
    d =  -0.87493628535349455 + 1j*0.48423805774828182;
    e =   0.90509035674203486 - 1j*0.42521929181608881;

    f = -(a + c + e + d + 1 + c2 + c/c1);
    b =  (a/c2 - c1/e) / (f' - d');
    z =  (c/c2 - c/e) / (c1/a - f/d);
    x = -(1 + c/c2 + d + z + b) / (b/z + b + 1);

    Y = [
        1  1     1      1    1       1      1       1    ;
        1  c1    a     -a    b      -1     -c1     -b    ;
        1  c     c2     a    d       c/c1   e       f    ;
        1  c/d   c      c/f  c/c2    c/a    c1      c/e  ;
        1  c/c2  d      z    b*x     x      b*x/z   b    ;
        1  c/f   e     -c/f  x      -x     -e      -1    ;
        1  c/e   c/c1  -z    z      -c/c1  -1      -c/e  ;
        1  c/a   f     -1    b*x/z  -c/a   -b*x/z  -f    ;
    ];

end
