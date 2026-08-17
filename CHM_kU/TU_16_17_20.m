function Y = TU_16_17_20(p)
% ------------------------------------------------------------------------------
% 2022-12-29 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 1-parameter 2-unitary CHM family.
% Looks like it cannot serve as self-dual.
%
% p1    d    #L   q : BH(16, q)
% ------------------------------------------------------------------------------
% r     17   20   --
% 1/5   17   20   20
% 1/3   17   12   12
% 1/7   17   20   28
% 1/9   17   20   36
% 3/16  17   16   16
% a/b   17   --   f(a,b) : f to be determined...
% 0     21   4    4
% 1/8   23   8    24
% ------------------------------------------------------------------------------

    a = exp(2j*pi*p);
    b = a^2;
    Y = [
        1   1   1       1     1     1     1     1     1     1     1     1       1     1     1     1    ;
        1  -1  -I*a     a    -a     I*a   1    -1     1    -1    -a     I*a     1    -1     a    -I*a  ;
        1   I   a      -a     a    -a    -I    -1    -I    -1     a    -a       I     1     a    -a    ;
        1   1  -1      -1    -1    -1     1     1     1     1    -1    -1       1     1    -1    -1    ;
        1  -1   1      -1    -1     1    -1     1     1/a   I/a  -I*a   a      -1/a  -I/a   I*a  -a    ;
        1   1   I*a     a    -a    -I*a  -1    -1    -I     I    -I*a  -a      -I     I     I*a   a    ;
        1  -I  -a      -a     a     a     I    -1    -1     I     I*a  -I*a     1    -I     I*a  -I*a  ;
        1  -1  -1       1     1    -1    -1     1    -1/a  -I/a  -I*a   a       1/a   I/a   I*a  -a    ;
        1  -I  -I*b    -b     I*a  -I*a   I*a  -I*a   I    -1     b     I*b    -I*a   I*a  -I*a   I*a  ;
        1  -1   I*a    -a    -I*a   a    -I    -I    -1     1    -a     I*a     I     I    -I*a   a    ;
        1   I  -a       a     I*a   I*a  -1    -I     I     1     a    -a      -1    -I    -I*a  -I*a  ;
        1  -I   I*b     b     I*a  -I*a  -I*a   I*a   I    -1    -b    -I*b     I*a  -I*a  -I*a   I*a  ;
        1   1   1       1    -I*a  -a     1/a  -I/a  -1/a   I/a   I*a   a      -1    -1    -1    -1    ;
        1  -1  -I*a     a    -I*a   a     I     I    -I    -I     I*a  -a      -1     1    -a     I*a  ;
        1   I   a      -a     I*a   I*a   1     I    -1    -I    -I*a  -I*a    -I    -1    -a     a    ;
        1   1  -1      -1    -I*a  -a    -1/a   I/a   1/a  -I/a   I*a   a      -1    -1     1     1    ;
    ];

    w = exp(2j*pi/12);
    DL = diag([1,1,1,1,1,1,-1,1,1,a,-a,-1,1,-1,-1,1]);
    DR = diag([1,1,1,1,1,-I,-I*b,-b,w,w^4,w^7,w^10,-I,-1,I,1]);
    Y = DL * Y * DR;

end

