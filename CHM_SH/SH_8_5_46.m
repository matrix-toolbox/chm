function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_5_46()
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 8 with generic d = 5 and #L = 46.
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 4;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.005;
end

function Y = pattern(p)
    Y = [
        1    1    1    1    1    1    1    1;
        1    1    1   -1   -1    1   -1   -1;
        1    1   -1    p(1)   -p(1)   -1    p(1)   -p(1);
        1   -1    p(1)    p(2)    p(3)   -p(1)   -p(2)   -p(3);
        1   -1   -p(1)    p(3)    p(4)    p(1)   -p(3)   -p(4);
        1    1   -1   -p(1)    p(1)   -1   -p(1)    p(1);
        1   -1    p(1)   -p(2)   -p(3)   -p(1)    p(2)    p(3);
        1   -1   -p(1)   -p(3)   -p(4)    p(1)    p(3)    p(4);
    ];
end
