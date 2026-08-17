function [iiMax, nP, ZTolerance, YPattern, muFactor] = DS_7()
% ------------------------------------------------------------------------------
% 2023-02-27 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% General pattern for doubly symmetric core of CHM of size N = 7.
%
% Recovers everything but C7A, C7B. :/
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 12;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.001;
end

function Y = pattern(p)
    Y = [
        1  1     1      1      1      1      1    ;
        1  p(1)  p(2)   p(3)   p(4)   p(5)   p(6) ;
        1  p(2)  p(7)   p(8)   p(9)   p(10)  p(5) ;
        1  p(3)  p(8)   p(11)  p(12)  p(9)   p(4) ;
        1  p(4)  p(9)   p(12)  p(11)  p(8)   p(3) ;
        1  p(5)  p(10)  p(9)   p(8)   p(7)   p(2) ;
        1  p(6)  p(5)   p(4)   p(3)   p(2)   p(1) ;
    ];
end
