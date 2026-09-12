function H=CBH9
% ------------------------------------------------------------------------------
% 2021-05-04 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric BH(9, 3);
% [1] https://math.uchicago.edu/~tghyde/Hyde,%20Kraisler%20--%20Circulant%20q-Butson%20Hadamard%20matrices.pdf
%
% >> H = CBH9;
% ------------------------------------------------------------------------------

    w = exp(2j*pi/3);

    H = [
        1 w^2 w 1 1 1 1 w w^2;
        w^2 1 w^2 w 1 1 1 1 w;
        w w^2 1 w^2 w 1 1 1 1;
        1 w w^2 1 w^2 w 1 1 1;
        1 1 w w^2 1 w^2 w 1 1;
        1 1 1 w w^2 1 w^2 w 1;
        1 1 1 1 w w^2 1 w^2 w;
        w 1 1 1 1 w w^2 1 w^2;
        w^2 w 1 1 1 1 w w^2 1;
    ];

end
