function [iiMax, nP, ZTolerance, YPattern, muFactor] = HH_10_6_106()
% ------------------------------------------------------------------------------
% 2023-01-23 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian CHM of size N = 10 with d = 6 and #L = 106.
% Original output from the Sinkhorn procedure with "unsorted" diagonal:
%
% 1   1       1     1          1       1          1        1           1       1        ;
% 1  -1       a     b^2/a      b^2     b          c        b           b^2/c  -b^2      ;
% 1   1/a     1    -1/a        c       a/c       -c       -a/c        -1      -1        ;
% 1   a/b^2  -a    -1          a       a*c/b^2/d  a*c/b^2 -a*c/b^2/d  -a/b^2  -a*c/b^2  ;
% 1   1/b^2   1/c   1/a       -1       1/b        c/b^2    1/b         a/b^2  -1/b^2    ;
% 1   1/b     c/a   d*b^2/c/a  b      -1          d        e           c/a    -1        ;
% 1   1/c    -1/c   b^2/c/a    b^2/c   1/d       -1       -1/d        -b^2/c  -b^2/c/a  ;
% 1   1/b    -c/a  -d*b^2/c/a  b       1/e       -d        1          -c/a    -1        ;
% 1   c/b^2  -1    -b^2/a      b^2/a   a/c       -c/b^2   -a/c         1      -1        ;
% 1  -1/b^2  -1    -b^2/c/a   -b^2    -1         -a*c/b^2 -1          -1       1        ;
% ------------------------------------------------------------------------------
% Call:
% >> Y = solver(nthargout([1:5], @HH_10_6_106)); summary(Y, 1e-8);
% ------------------------------------------------------------------------------

    iiMax = 10000;
    nP = 5;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.01;
end

function Y = pattern(p)
    Y = [
        1   1           1                       1              1                  1             1                      1             1                       1                ;
        1   1          -p(1)/p(3)              -1             -1                  1/p(1)       -1/p(1)                 p(3)          p(1)/p(3)              -p(3)             ;
        1  -p(3)/p(1)   1                      -p(3)/p(1)     -1                  1/p(2)       -p(4)*p(2)^2/p(3)/p(1)  p(2)          1/p(5)                 -p(4)             ;
        1  -1          -p(1)/p(3)               1             -1                  p(3)/p(2)^2  -p(2)^2/p(1)            p(2)^2/p(1)   p(1)/p(3)              -p(3)/p(2)^2      ;
        1  -1          -1                      -1              1                 -1/p(2)^2     -p(2)^2/p(3)/p(1)      -p(2)^2       -1                      -p(1)*p(3)/p(2)^2 ;
        1   p(1)        p(2)                    p(2)^2/p(3)   -p(2)^2            -1             p(2)^2/p(1)            p(2)^2        p(2)                    p(3)             ;
        1  -p(1)       -p(1)*p(3)/p(2)^2/p(4)  -p(1)/p(2)^2   -p(1)*p(3)/p(2)^2   p(1)/p(2)^2  -1                      p(1)          p(1)*p(3)/p(2)^2/p(4)   p(1)*p(3)/p(2)^2 ;
        1   1/p(3)      1/p(2)                  p(1)/p(2)^2   -1/p(2)^2           1/p(2)^2      1/p(1)                -1             1/p(2)                  p(3)/p(2)^2      ;
        1   p(3)/p(1)   p(5)                    p(3)/p(1)     -1                  1/p(2)        p(4)*p(2)^2/p(3)/p(1)  p(2)         -1                       p(4)             ;
        1  -1/p(3)     -1/p(4)                 -p(2)^2/p(3)   -p(2)^2/p(3)/p(1)   1/p(3)        p(2)^2/p(3)/p(1)       p(2)^2/p(3)   1/p(4)                 -1                ;
    ];
end


