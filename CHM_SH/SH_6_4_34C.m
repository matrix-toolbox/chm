function Y = SH_6_4_34C
% ------------------------------------------------------------------------------
% 2026-09-10
%
% symmetric CHM(6) with d = 4 and L = 34
% this is a companion ("C") solution to SH_6_4_34.m
% -- same order, defect and Haagerup count, but verified monomially INEQUIVALENT
%
% c is fixed to a primitive 5th root of unity
% b solves the quadratic found by (Wolfram's) -- both roots provide monomially equivalent matrices
%
% symmetric-pattern consistency condition ==> Y_22 = a*c^2/b = -1 - a - b - 2*c
% ------------------------------------------------------------------------------

    c = exp(2j*pi/5);

    b = ( -c*(1+c)^2 + sqrt(2)*sqrt(c^2*(1+2*c+2*c^3+c^4)) ) / (c*(2+c) - 1);
    a = -(b + b^2 + 2*b*c) / (b + c^2);

    Y = [
        1        1        1          1     1    1   ;
        1  a*c^2/b        a          b     c    c   ;
        1        a       -a        b/c  -b/c   -1   ;
        1        b      b/c  b^2/c^2/a   b/c  b/a   ;
        1        c     -b/c        b/c    -1   -c   ;
        1        c       -1        b/a    -c  -b/a  ;
    ];

end
