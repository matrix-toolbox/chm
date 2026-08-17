function Y = Y_8_0_10
% ------------------------------------------------------------------------------
% 2022-08-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% An isolated matrix with #L = 10 found by the Sinkhorn algorithm.
% ------------------------------------------------------------------------------

    a = 1/3 - 2j*sqrt(2)/3;

    Y = [
        1      1      1      1      1      1      1      1     ;
        1      1     -1     -1      a     -a      a'    -a'    ;
        1      1     -a'     a'    -a      a     -1     -1     ;
        1     -a      a'     1      a     -1     -a'    -1     ;
        1     -a     -1     -a'    -1      a      1      a'    ;
        1      a^2    1     -a      a^2    a^2   -a      1     ;
        1      a^2    a     -1     -a^2   -a      a     -a     ;
        1      a^2   -a      a     -a     -a^2   -1      a     ;
    ];

end
