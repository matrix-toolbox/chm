function [iiMax, nP, ZTolerance, YPattern, muFactor] = HH_12_9_272()
% ------------------------------------------------------------------------------
% 2023-01-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian CHM of size N = 12 with generic d = 9 and #L = 272.
% Found as a solution of "HH_12_5_1902.m" and defined as a separate pattern.
% ------------------------------------------------------------------------------

    iiMax = 50000;
    nP = 18;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.0004;
end

function Y = pattern(p)
    Y = [
        1   1      1       1       1       1      1       1        1       1       1       1      ;
        1   1      p(1)    p(2)    p(1)    p(2)   p(3)    p(4)     p(5)    p(6)    p(7)    p(8)   ;
        1   p(1)'  1       p(9)    p(1)'   p(2)   p(10)   p(11)    p(12)   p(13)   p(14)   p(15)  ;
        1   p(2)'  p(9)'   1       p(1)    p(1)   p(15)'  p(12)'   p(14)'  p(10)'  p(11)'  p(13)' ;
        1   p(1)'  p(1)    p(1)'   1       p(1)   p(16)   p(16)'   p(17)   p(18)   p(18)'  p(17)' ;
        1   p(2)'  p(2)'   p(1)'   p(1)'   1      p(5)'   p(8)'    p(6)'   p(4)'   p(3)'   p(7)'  ;
        1   p(3)'  p(10)'  p(15)   p(16)'  p(5)  -1      -p(16)'  -p(5)   -p(10)' -p(3)'  -p(15)  ;
        1   p(4)'  p(11)'  p(12)   p(16)   p(8)  -p(16)  -1       -p(12)  -p(4)'  -p(11)' -p(8)   ;
        1   p(5)'  p(12)'  p(14)   p(17)'  p(6)  -p(5)'  -p(12)'  -1      -p(6)   -p(14)  -p(17)' ;
        1   p(6)'  p(13)'  p(10)   p(18)'  p(4)  -p(10)  -p(4)    -p(6)'  -1      -p(18)' -p(13)' ;
        1   p(7)'  p(14)'  p(11)   p(18)   p(3)  -p(3)   -p(11)   -p(14)' -p(18)  -1      -p(7)'  ;
        1   p(8)'  p(15)'  p(13)   p(17)   p(7)  -p(15)' -p(8)'   -p(17)  -p(13)  -p(7)   -1      ;
    ];
end
