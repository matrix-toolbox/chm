function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_7_general()
% ------------------------------------------------------------------------------
% 2023-01-31 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% General pattern for symmetric CHM of size N = 7.
% It should reproduce all other cases.
% When on a right path, it returns result immediately,
% but usually gets stuck "in the middle" and requires a nudge.
%
% Found solutions:
% (d, #L) in {
%              (0, 5)
%              (0, 7) : BH(7, 7)
%              (3, 6) : BH(7, 6)
%            }
%
% In particular, never observed the maximal value of #L = 463.
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 21;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.002;
end

function Y = pattern(p)
    Y = [
        1  1     1      1      1      1      1     ;
        1  p(1)  p(2)   p(3)   p(4)   p(5)   p(6)  ;
        1  p(2)  p(7)   p(8)   p(9)   p(10)  p(11) ;
        1  p(3)  p(8)   p(12)  p(13)  p(14)  p(15) ;
        1  p(4)  p(9)   p(13)  p(16)  p(17)  p(18) ;
        1  p(5)  p(10)  p(14)  p(17)  p(19)  p(20) ;
        1  p(6)  p(11)  p(15)  p(18)  p(20)  p(21) ;
    ];
end
