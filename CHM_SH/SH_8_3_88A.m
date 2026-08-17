function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_3_88A()
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 8 with generic d = 3 and #L = 88.
% Compare with "SH_8_3_50.m" and "SH_8_3_88B.m".
% ------------------------------------------------------------------------------
% It contains 1-parameter family with d = 11 and #L = 6, see "SH_8_11_6A.m".
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 13;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.005;
end

function Y = pattern(p)
    Y = [
        1    1       1       1       1       1       1            1     ;
        1    p(1)    p(2)   -1       p(3)   -p(3)   -p(1)        -p(2)  ;
        1    p(2)    p(4)    p(5)   -1      -p(4)   -p(5)        -p(2)  ;
        1   -1       p(5)    p(6)    p(7)   -p(6)   -p(5)        -p(7)  ;
        1    p(3)   -1       p(7)    p(8)   -p(3)   -p(8)        -p(7)  ;
        1   -p(3)   -p(4)   -p(6)   -p(3)    p(9)    p(10)        p(11) ;
        1   -p(1)   -p(5)   -p(5)   -p(8)    p(10)   p(10)/p(11)  p(12) ;
        1   -p(2)   -p(2)   -p(7)   -p(7)    p(11)   p(12)        p(13) ;
    ];
end
