function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_9_0_76()
% ------------------------------------------------------------------------------
% 2023-01-15 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 9 with generic d = 0 and #L = 76.
% There might be some additional dependencies between entries of Y.
% ------------------------------------------------------------------------------

    iiMax = 10000;
    nP = 7;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.01;
end

function Y = pattern(p)
    x = -(1+2*p(1)+2*p(2)+p(3)+2*p(5)/p(4));
    Y = [
        1     1           1            1                1                 1         1          1                   1                   ;
        1     x           p(1)         p(2)             p(2)              p(3)      p(1)       p(5)/p(4)           p(5)/p(4)           ;
        1     p(1)        p(5)        -1                p(6)             -p(5)      p(4)       p(5)/p(4)           p(5)/p(6)           ;
        1     p(2)       -1            p(7)             p(4)*p(7)/p(5)   -p(7)      p(6)       p(5)/p(4)           p(7)/p(6)           ;
        1     p(2)        p(6)         p(4)*p(7)/p(5)   p(7)             -p(7)     -1          p(7)/p(6)           p(5)/p(4)           ;
        1     p(3)       -p(5)        -p(7)            -p(7)              p(3)^2   -p(5)       p(3)                p(3)                ;
        1     p(1)        p(4)         p(6)            -1                -p(5)      p(5)       p(5)/p(6)           p(5)/p(4)           ;
        1     p(5)/p(4)   p(5)/p(4)    p(5)/p(4)        p(7)/p(6)         p(3)      p(5)/p(6)  p(5)^2/p(4)^2       p(5)^2/p(4)^2/p(3)  ;
        1     p(5)/p(4)   p(5)/p(6)    p(7)/p(6)        p(5)/p(4)         p(3)      p(5)/p(4)  p(5)^2/p(4)^2/p(3)  p(5)^2/p(4)^2       ;
    ];
end
