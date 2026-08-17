function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_3_142()
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 8 with generic d = 3 and #L = 142.
%
% It also contains solutions for several families:
% * SH_8_3_46(p1, p2)
% * SH_8_5_26(p1, p2)
% * SH_8_9_10(p1)
% * SH_8_13_6A(p1)
% * SH_8_13_6B(p1)
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 9;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.005;
end

function Y = pattern(p)
    x1 = -(1+p(1)*p(4)/p(9)+p(1)/p(3)+p(1)+p(2)+p(9)/p(4)+p(3));
    x2 = -(1+p(1)/p(3)-p(1)*p(4)/p(9)+p(5)+p(6)-p(7)/p(4)-p(2));
    x3 = -(1+p(1)+p(1)*p(4)^2/p(9)+p(5)+p(7)+p(8)+p(9));

    Y = [
        1   1                1                  1               1                 1          1           1         ;
        1   x1               p(1)*p(4)/p(9)     p(1)/p(3)       p(1)              p(2)       p(9)/p(4)   p(3)      ;
        1   p(1)*p(4)/p(9)  -p(1)*p(4)^2/p(9)  -p(1)*p(4)/p(9)  p(1)*p(4)^2/p(9)  p(4)      -1          -p(4)      ;
        1   p(1)/p(3)       -p(1)*p(4)/p(9)     x2              p(5)              p(6)      -p(7)/p(4)  -p(2)      ;
        1   p(1)             p(1)*p(4)^2/p(9)   p(5)            x3                p(7)       p(8)        p(9)      ;
        1   p(2)             p(4)               p(6)            p(7)              p(7)/p(2)  p(7)/p(4)   p(9)/p(3) ;
        1   p(9)/p(4)       -1                 -p(7)/p(4)       p(8)              p(7)/p(4) -p(8)       -p(9)/p(4) ;
        1   p(3)            -p(4)              -p(2)            p(9)              p(9)/p(3) -p(9)/p(4)  -p(9)/p(2) ;
    ];
end
