function [iiMax, nP, ZTolerance, YPattern, muFactor] = HH_12_9_314()
% ------------------------------------------------------------------------------
% 2023-01-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian CHM of size N = 12 with generic d = 9 and #L = 314.
% Found as a solution of "HH_12_5_1902.m" and defined as a separate pattern.
% ------------------------------------------------------------------------------

    iiMax = 50000;
    nP = 8;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.0004;
end

function Y = pattern(p)
    Y = [
        1   1               1           1               1                1                1                1               1      1               1                1               ;
        1   1              -1           1j*p(2)        -p(1)*p(4)/p(5)   p(1)             p(2)             p(3)           -p(3)   p(4)            1j               p(5)            ;
        1  -1               1          -1               p(4)/p(5)        p(5)/p(4)       -1                1               1     -p(4)/p(5)      -1               -p(5)/p(4)       ;
        1  -1j/p(2)        -1           1               p(6)            -p(6)*p(5)/p(4)   1j              -p(3)            p(3)   p(7)           -p(2)'            p(7)*p(5)/p(4)  ;
        1  -p(5)/p(1)/p(4)  p(5)/p(4)   p(6)'           1                p(5)/p(4)       -p(6)'           -p(8)*p(5)/p(4)  p(8)  -p(8)'           p(5)/p(1)/p(4)   p(5)/p(4)/p(8)  ;
        1   p(1)'           p(4)/p(5)  -p(4)/p(5)/p(6)  p(4)/p(5)        1                p(4)/p(5)/p(6)   p(4)*p(8)/p(5) -p(8)  -p(4)/p(5)/p(8) -p(1)'            p(8)'           ;
        1   p(2)'          -1          -1j             -p(6)             p(6)*p(5)/p(4)  -1               -p(3)            p(3)  -p(7)            1j/p(2)         -p(7)*p(5)/p(4)  ;
        1   p(3)'           1          -p(3)'          -p(4)/p(5)/p(8)   p(5)/p(4)/p(8)  -p(3)'           -1              -1      p(4)/p(5)/p(8)  p(3)'           -p(5)/p(4)/p(8)  ;
        1  -p(3)'           1           p(3)'           p(8)'           -p(8)'            p(3)'           -1              -1      p(8)'          -p(3)'           -p(8)'           ;
        1   p(4)'          -p(5)/p(4)   p(7)'          -p(8)            -p(8)*p(5)/p(4)  -p(7)'            p(8)*p(5)/p(4)  p(8)  -1              -p(4)'            p(5)/p(4)       ;
        1  -1j             -1          -p(2)            p(1)*p(4)/p(5)  -p(1)            -1j*p(2)          p(3)           -p(3)  -p(4)           -1               -p(5)            ;
        1   p(5)'          -p(4)/p(5)   p(4)/p(5)/p(7)  p(4)*p(8)/p(5)   p(8)            -p(4)/p(5)/p(7)  -p(4)*p(8)/p(5) -p(8)   p(4)/p(5)      -p(5)'           -1               ;
    ];
end
