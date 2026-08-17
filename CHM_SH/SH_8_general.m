function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_general()
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% General pattern for symmetric CHM of order N = 8
% with generic d = 5 and max. possible #L = 813.
% ------------------------------------------------------------------------------
% Convergence is quite slow and works only for particular initial conditions.
% It should reproduce all other cases.
% Generic values of (d, #L) are (5, 813) as expected.
% It generates also:
% (d, #L) in { (0, 10), (3, 76), (3, 88), (3, 142),  (5, 813), ... }
% and for p(1) = -1 ==> (5, 46), (3, 76), ...
% ------------------------------------------------------------------------------


    iiMax = 40000;
    nP = 21;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.001;
end

function Y = pattern(p)
    x1 = -(1+p(1)+p(2)+p(3)+p(4)+p(5)+p(6));
    x2 = -(1+p(1)+p(7)+p(8)+p(9)+p(10)+p(11));
    x3 = -(1+p(2)+p(7)+p(12)+p(13)+p(14)+p(15));
    x4 = -(1+p(3)+p(8)+p(12)+p(16)+p(17)+p(18));
    x5 = -(1+p(4)+p(9)+p(13)+p(16)+p(19)+p(20));
    x6 = -(1+p(5)+p(10)+p(14)+p(17)+p(19)+p(21));
    x7 = -(1+p(6)+p(11)+p(15)+p(18)+p(20)+p(21));

    Y = [
        1   1      1      1       1      1       1       1     ;
        1   x1     p(1)   p(2)    p(3)   p(4)    p(5)    p(6)  ;
        1   p(1)   x2     p(7)    p(8)   p(9)    p(10)   p(11) ;
        1   p(2)   p(7)   x3      p(12)  p(13)   p(14)   p(15) ;
        1   p(3)   p(8)   p(12)   x4     p(16)   p(17)   p(18) ;
        1   p(4)   p(9)   p(13)   p(16)  x5      p(19)   p(20) ;
        1   p(5)   p(10)  p(14)   p(17)  p(19)   x6      p(21) ;
        1   p(6)   p(11)  p(15)   p(18)  p(20)   p(21)   x7    ;
    ];
end
