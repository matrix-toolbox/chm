function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_3_50()
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of size N = 8 with generic d = 3 and #L = 50.
% Constructed by a mistake! :)
% This should be "SH_8_3_88X.m", but some parameters were accidentally mixed...
% Nevertheless, RW finds a solution.
%
% Special cases:
% * p(1) = -1            ==> contains BH(8, 2) with p2+ = -+-+++-+-++ and d = 21
% * p(1) = -1            ==> contains non-BH with d = 3 and #L = 18
% * p(1) = exp(2j*pi/11) ==> generic non-BH solution with (d, #L) = (3, 50)
% * p(1) = -1j;          ==> non-BH with d = 3, #L = 48
%                            p1 =  -  1i
%                            p2 = -0.973797648849187 - 0.227416224346013i
%                            p3 = +0.964937884619837 + 0.262478339726147i
%                            p4 =  -  1i
%                            p5 = +0.964937884619819 + 0.262478339726213i
%                            p6 = -0.506550587787758 + 0.862210242349206i
%                            p7 = +0.724420484698463 + 0.689358369318342i
%                            p8 = -0.506550587787725 + 0.862210242349225i
%                            p9 = +0.862210242349173 + 0.506550587787815i
%                           p10 = -0.475043320547812 + 0.879962410335185i
%                           p11 = +0.036158259650669 - 0.999346076321429i
%                           p12 = -0.548667687205856 - 0.836040530725740i
% * A8 can also be recovered from this pattern...
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 12;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.005;
end

function Y = pattern(p)
    Y = [
        1       1       1       1      1       1       1           1     ;
        1       p(1)    p(2)   -1      p(3)   -p(3)   -p(1)       -p(2)  ;
        1       p(2)    p(4)    p(5)  -1      -p(4)   -p(5)       -p(2)  ;
        1      -1       p(5)    p(6)   p(7)   -p(6)   -p(5)       -p(7)  ;
        1       p(3)   -1       p(7)   p(8)   -p(3)   -p(8)       -p(7)  ;
        1      -p(3)   -p(4)   -p(6)  -p(3)    p(11)   p(9)        p(10) ;
        1      -p(1)   -p(5)   -p(5)  -p(8)    p(9)    p(9)/p(10)  p(10) ;
        1      -p(2)   -p(2)   -p(7)  -p(7)    p(10)   p(10)       p(12) ;
    ];
end
