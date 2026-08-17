function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_3_88B()
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 8 with generic d = 3 and #L = 88.
% Compare with "SH_8_3_50.m" and "SH_8_3_88A.m".
% ------------------------------------------------------------------------------
% It contains 1-parameter family with d = 11 and #L = 6, see "SH_8_11_6B.m".
% ..............................................................................
% Other possible solutions can be obtained with fixed p(1):
% p(1) = -1 ==> d = 3, #L = 66, non-BH (unresolved case)
% p(1) = -I ==> non-BH, d = 3, #L = 86 (unresolved case)
%               or BH(8, 4), d = 11, #L = 4, see "BH_8_11_4.m";
%
%               0  0  0  0  0  0  0  0
%               0  3  0  1  2  2  0  2
%               0  0  3  2  2  1  2  0
%               0  1  2  3  2  0  0  2
%               0  2  2  2  0  2  0  0
%               0  2  1  0  2  3  2  0
%               0  0  2  0  0  2  2  2
%               0  2  0  2  0  0  2  2
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 9;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.005;
end

function Y = pattern(p)
    x = -(1-p(1)+p(5)+2*p(4)+p(3)+p(9));
    y = -(1+2*p(3)/p(4)+p(7)+p(3)-p(8)+p(3)/p(9));
    Y = [
        1    1          1          1       1       1         1           1         ;
        1    p(1)       p(2)      -p(1)   -1       p(3)/p(4) -p(3)/p(4) -p(2)      ;
        1    p(2)       p(5)*p(7)  p(5)    p(6)    p(7)       p(6)       p(2)      ;
        1   -p(1)       p(5)       x       p(4)    p(3)       p(9)       p(4)      ;
        1   -1          p(6)       p(4)    p(8)   -p(8)      -p(6)      -p(4)      ;
        1    p(3)/p(4)  p(7)       p(3)   -p(8)    y          p(3)/p(4)  p(3)/p(9) ;
        1   -p(3)/p(4)  p(6)       p(9)   -p(6)    p(3)/p(4) -p(9)      -1         ;
        1   -p(2)       p(2)       p(4)   -p(4)    p(3)/p(9) -1         -p(3)/p(9) ;
    ];
end

