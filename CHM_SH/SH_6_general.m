function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_6_general()
% ------------------------------------------------------------------------------
% 2023-02-02 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% The most general pattern for a symmetric CHM of size N = 6
% with d = 4 and max. possible #L = 241.
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 15;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.001;
end

function Y = pattern(p)
    Y = [
        1  1     1     1      1      1     ;
        1  p(1)  p(2)  p(3)   p(4)   p(5)  ;
        1  p(2)  p(6)  p(7)   p(8)   p(9)  ;
        1  p(3)  p(7)  p(10)  p(11)  p(12) ;
        1  p(4)  p(8)  p(11)  p(13)  p(14) ;
        1  p(5)  p(9)  p(12)  p(14)  p(15) ;
    ];
end
