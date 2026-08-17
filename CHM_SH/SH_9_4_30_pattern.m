function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_9_4_30_pattern()
% ------------------------------------------------------------------------------
% 2023-01-15 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 9 with generic d = 4 and #L = 30.
%
% Parameter p(2) might be "anything":
%    p(2) = +1 ==> BH(9, 6) with d = 10, #L = 6 (unresolved case)
%    p(2) = +I ==> BH(9, 12) with d = 4, #L = 12 (unresolved case)
%
% Indeed, it contains 1-parameter family, see "SH_9_4_30.m".
% ------------------------------------------------------------------------------

    iiMax = 10000;
    nP = 3;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.0002;
end

function Y = pattern(p)
    w = exp(2j*pi/3);
    Y = [
        1    1       1    1       1       1       1       1       1     ;
        1    w      -w'   p(1)    p(2)   -p(1)    w'     -p(2)    w'    ;
        1   -w'      w    w'      w'      w'     -w'      w'     -w'    ;
        1    p(1)    w'   w       w'      w       p(3)    1       p(2)' ;
        1    p(2)    w'   w'      w       1       p(1)'   w       p(3)' ;
        1   -p(1)    w'   w       1       w      -p(3)    w'     -p(2)' ;
        1    w'     -w'   p(3)    p(1)'  -p(3)    w      -p(1)'   w'    ;
        1   -p(2)    w'   1       w       w'     -p(1)'   w      -p(3)' ;
        1    w'     -w'   p(2)'   p(3)'  -p(2)'   w'     -p(3)'   w     ;
    ];

end
