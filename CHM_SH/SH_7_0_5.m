function H = SH_7_0_5
% ------------------------------------------------------------------------------
% 2023-01-16 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of size N = 7 with d = 0 and #L = 5.
% ------------------------------------------------------------------------------

    a = -3/4 - 1j*sqrt(7)/4;
    b = 1/a^2;
    H = [
        1       1       1       1       1       1       1   ;
        1       a       a       1/a     1       1/a     1   ;
        1       a       1       1       a       1/a     1/a ;
        1       1/a     1       1/a     1/a     b       b   ;
        1       1       a       1/a     a       1       1/a ;
        1       1/a     1/a     b       1       b       1/a ;
        1       1       1/a     b       1/a     1/a     b   ;
    ];

end
