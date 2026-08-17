function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_10_0_134()
% ------------------------------------------------------------------------------
% 2023-01-21 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 10 with generic d = 0 and #L = 134.
%
% It contains solutions:
% s1) non-BH with d = 0 and #L = 134 (generic).
% s2) BH(10, 10) with d = 8 and #L = 10; see: "BH_10_8_10.m".
% s3) non-BH with d = 8, #L = 58; see "SH_10_8_58.m".
% ----------------------------------------------------------------------

    iiMax = 40000;
    nP = 6;
    ZTolerance = 6e-13;
    YPattern = @pattern;
    muFactor = 0.003;
end

function Y = pattern(p)
    x = p(2)*p(5)/p(1)/p(6);
    y = -p(5)*p(3)^2*x^2/p(2)^3;
    Y = [
        1   1             1             1                1             1             1       1               1                1            ;
        1   p(2)^2/p(1)   p(1)          p(2)             p(3)          p(4)          x       p(2)^2/x        p(2)^2          -p(2)         ;
        1   p(1)         -y*p(5)*x     -p(2)*p(5)*y     -p(1)          p(5)          y       p(2)/x         -p(2)*y          -p(2)*p(5)    ;
        1   p(2)         -p(2)*p(5)*y  -p(5)*y           p(2)*y/x      p(2)*p(5)    -y       p(1)/x          p(2)*p(5)/p(6)  -p(5)*x       ;
        1   p(3)         -p(1)          p(2)*y/x         p(3)/x/y      p(6)         -1/x     p(2)^2/x^2     -p(2)^2/x         p(2)/x/y     ;
        1   p(4)          p(5)          p(2)*p(5)        p(6)          p(5)*x/y     -1/y    -p(2)/x         -p(2)/y          -p(2)*p(5)/y  ;
        1   x             y            -y               -1/x          -1/y          -1       1/x            -x                1/y          ;
        1   p(2)^2/x      p(2)/x        p(1)/x           p(2)^2/x^2   -p(2)/x        1/x     p(3)/x          p(2)^2/x/p(1)   -p(1)/y       ;
        1   p(2)^2       -p(2)*y        p(2)*p(5)/p(6)  -p(2)^2/x     -p(2)/y       -x       p(2)^2/x/p(1)  -p(3)*x           p(1)/y       ;
        1  -p(2)         -p(2)*p(5)    -p(5)*x           p(2)/x/y     -p(2)*p(5)/y   1/y    -p(1)/y          p(1)/y           p(5)/y       ;
    ];
end
