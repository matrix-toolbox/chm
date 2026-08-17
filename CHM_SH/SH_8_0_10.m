function Y = SH_8_0_10
% ------------------------------------------------------------------------------
% 2023-01-16 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of size N = 8 with d = 0 and #L = 10.
% ------------------------------------------------------------------------------

    a = 1/3 - 1j*sqrt(8)/3;
    b = -a^2;
    Y = [
        1    1    1    1    1    1    1    1  ;
        1    a   -1    b   -b    a   -a   -a  ;
        1   -1    a'  -1    1    a   -a   -a' ;
        1    b   -1   -a   -a   -b    a    a  ;
        1   -b    1   -a   -b   -a   -b    1  ;
        1    a    a   -b   -a   -a    b   -1  ;
        1   -a   -a    a   -b    b    a   -1  ;
        1   -a   -a'   a    1   -1   -1    a' ;
    ];

end
