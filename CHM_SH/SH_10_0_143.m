function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_10_0_143()
% ------------------------------------------------------------------------------
% 2023-01-21 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 10 with generic d = 0 and #L = 143.
%
% It is a special solution of "SH_10_0_349A.m".
% ----------------------------------------------------------------------

    iiMax = 30000;
    nP = 16;
    ZTolerance = 6e-13;
    YPattern = @pattern;
    muFactor = 0.003;
end

function Y = pattern(p)
    x = -(1+4*p(5)+2*p(10)+2*p(13));
    y = -(2+4*p(8)+p(5)+2*p(15));
    Y = [
        1  1     1       1       1      1       1       1       1       1   ;
        1  p(1)  p(2)    p(3)    p(4)   p(5)    p(6)    p(7)    p(8)    p(9)   ;
        1  p(2)  p(2)    p(3)    1      p(10)   p(11)   p(3)    p(8)    p(12)   ;
        1  p(3)  p(3)    p(2)    1      p(10)   p(12)   p(2)    p(8)    p(11)   ;
        1  p(4)  1       1       p(5)'  p(5)    p(4)'   p(4)    1       p(4)'  ;
        1  p(5)  p(10)   p(10)   p(5)   x       p(13)   p(5)    p(5)    p(13)   ;
        1  p(6)  p(11)   p(12)   p(4)'  p(13)   p(14)   p(9)    p(15)   p(16)   ;
        1  p(7)  p(3)    p(2)    p(4)   p(5)    p(9)    p(1)    p(8)    p(6)   ;
        1  p(8)  p(8)    p(8)    1      p(5)    p(15)   p(8)    y       p(15)   ;
        1  p(9)  p(12)   p(11)   p(4)'  p(13)   p(16)   p(6)    p(15)   p(14)   ;
    ];
end
