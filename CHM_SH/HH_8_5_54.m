function [iiMax, nP, ZTolerance, YPattern, muFactor] = HH_8_5_54()
% ------------------------------------------------------------------------------
% 2023-01-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian CHM of size N = 8 with d = 5 and #L = 54.
%
% Call:
% >> Y = solver(nthargout([1:5], @HH_8_5_54)); summary(Y, 1e-8);
% ------------------------------------------------------------------------------

    iiMax = 10000;
    nP = 3; % number of phases
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.01;
end

function Y = pattern(p)
    Y = [
        1  1          1               1           1                1           1      1                   ;
        1  1          p(1)            p(1)        p(2)             p(3)        p(1)/p(3)   p(1)/p(2)      ;
        1  p(1)'      1               p(1)        p(2)*p(3)/p(1)   p(3)/p(2)   p(2)/p(3)   p(1)/p(2)/p(3) ;
        1  p(1)'      p(1)'           1           p(3)/p(1)        p(2)'       p(2)/p(1)   p(3)'          ;
        1  p(2)'      p(1)/p(2)/p(3)  p(1)/p(3)  -1               -p(2)'      -p(1)/p(3)  -p(1)/p(2)/p(3) ;
        1  p(3)'      p(2)/p(3)       p(2)       -p(2)            -1          -p(2)/p(3)  -p(3)'          ;
        1  p(3)/p(1)  p(3)/p(2)       p(1)/p(2)  -p(3)/p(1)       -p(3)/p(2)  -1          -p(1)/p(2)      ;
        1  p(2)/p(1)  p(2)*p(3)/p(1)  p(3)       -p(2)*p(3)/p(1)  -p(3)       -p(2)/p(1)  -1              ;
    ];
end
