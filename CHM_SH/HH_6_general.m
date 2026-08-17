function [iiMax, nP, ZTolerance, YPattern, muFactor] = HH_6_general
% ------------------------------------------------------------------------------
% 2023-02-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Generic Hermitian case for N = 6, returns Y = Y^* : d(Y) = 4, #L(Y) = 34.
%
% Call:
% >> Y = solver(nthargout([1:5], @HH_6_general)); summary(Y, 1e-8);
% ------------------------------------------------------------------------------

    iiMax = 10000;
    nP = 15; % number of phases
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.006;
end

function Y = pattern(p)
    Y = [
        1    1      1      1      1       1      ;
        1    p(11)  p(1)   p(2)   p(3)    p(4)   ;
        1    p(1)'  p(12)  p(5)   p(6)    p(7)   ;
        1    p(2)'  p(5)'  p(13)  p(8)    p(9)   ;
        1    p(3)'  p(6)'  p(8)'  p(14)   p(10)  ;
        1    p(4)'  p(7)'  p(9)'  p(10)'  p(15)  ;
    ];
end
