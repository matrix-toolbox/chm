function H = C6
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% C-6-R CHM of order N = 6 with defect = 3;
% No affine family stems from C6, however it belongs to --> B6.m;
% [1] J. Symb. Comput. 12, 329--336, (1991)
% [2] https://arxiv.org/abs/math/0609076
%
% >> H = C6;
% ------------------------------------------------------------------------------

    d = (0.5 * (1 - sqrt(3))) + i * sqrt(0.5 * sqrt(3));

    H = [
        1       1       1       1       1     1;
        1      -1      -d    -d^2     d^2     d;
        1 -d^(-1)       1     d^2    -d^3   d^2;
        1 -d^(-2)  d^(-2)      -1     d^2  -d^2;
        1  d^(-2) -d^(-3)  d^(-2)       1    -d;
        1  d^(-1)  d^(-2) -d^(-2) -d^(-1)    -1;
    ];

end
