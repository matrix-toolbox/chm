function H = CBH8
% ------------------------------------------------------------------------------
% 2021-05-04 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric BH(8, 4);
% [1] https://math.uchicago.edu/~tghyde/Hyde,%20Kraisler%20--%20Circulant%20q-Butson%20Hadamard%20matrices.pdf
%
% >> H = CBH8;
% ------------------------------------------------------------------------------

    H = [
        1 -1  i  1  1  1  i -1;
       -1  1 -1  i  1  1  1  i;
        i -1  1 -1  i  1  1  1;
        1  i -1  1 -1  i  1  1;
        1  1  i -1  1 -1  i  1;
        1  1  1  i -1  1 -1  i;
        i  1  1  1  i -1  1 -1;
       -1  i  1  1  1  i -1  1;
    ];

end
