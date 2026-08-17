function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_8_0_70
% ------------------------------------------------------------------------------
% 2023-01-16 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of size N = 8 with d = 0 and #L = 70.
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 3;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.001;
end

function Y = pattern(p)
    g =  p(1)^2;
    k = -p(2)^2;
    d = -p(1)/p(2);
    b =  p(3)*p(1);
    f =  d/p(3);
    e =  p(1)^2/f;
    j = -f*p(2);
    m = -p(2)/b;
    o =  m*f;
    q =  g/d;

    Y = [
        1   1      1        1    1         1     1      1    ;
        1   p(3)   b       -1    p(3)'     d     d'     b'   ;
        1   b     -p(1)*d   e    f         g     p(2)   j    ;
        1  -1      e        k    m        -e    -k     -m    ;
        1   p(3)'  f        m    1/p(3)/b  p(1) -p(3)'  o    ;
        1   d      g       -e    p(1)     -p(1)  q     -f    ;
        1   d'     p(2)    -k   -p(3)'     q    -q/f   -p(2) ;
        1   b'     j       -m    o        -f    -p(2)   j/e  ;
    ];
end
