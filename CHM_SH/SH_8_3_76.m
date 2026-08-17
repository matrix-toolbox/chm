function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_3_76()
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of size N = 8 with d = 3 and #L = 76.
% ------------------------------------------------------------------------------

    iiMax = 10000;
    nP = 8;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.01;
end

function Y = pattern(p)
    Y = [
        1     1           1           1        1            1        1            1         ;
        1     p(1)       -1           p(2)     p(3)/p(4)   -p(1)    -p(3)/p(4)   -p(2)      ;
        1    -1           p(5)        p(4)     p(3)/p(2)   -p(5)    -p(4)        -p(3)/p(2) ;
        1     p(2)        p(4)        p(6)    -1           -p(6)    -p(2)        -p(4)      ;
        1     p(3)/p(4)   p(3)/p(2)  -1        p(7)        -p(7)    -p(3)/p(2)   -p(3)/p(4) ;
        1    -p(1)       -p(5)       -p(6)    -p(7)         p(8)^2   p(8)         p(8)      ;
        1    -p(3)/p(4)  -p(4)       -p(2)    -p(3)/p(2)    p(8)     p(3)         p(3)/p(8) ;
        1    -p(2)       -p(3)/p(2)  -p(4)    -p(3)/p(4)    p(8)     p(3)/p(8)    p(3)      ;
    ];
end
