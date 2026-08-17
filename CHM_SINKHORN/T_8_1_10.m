function Y=T_8_1_10(p1)
% ------------------------------------------------------------------------------
% 2022-04-09 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A Symmetric 1-parameter family obtained from a random seed applied to the Sinkhorn algorithm.
% Generic values: d(Y) = 3, #L(Y) = 10.
% It stems from BH(8, 2).
%
% [1] https://arxiv.org/abs/2204.11727
% ------------------------------------------------------------------------------

    y = exp(2j*pi*p1);

    Y = [
        1    1      1    1     1     1     1     1  ;
        1    y^2   -y^2 -y^2   y^2  -1    -1     1  ;
        1   -y^2   -y    y^2  -y     y    -1     y  ;
        1   -y^2    y^2  y     y    -1    -y    -y  ;
        1    y^2   -y    y    -y^2  -y     y    -1  ;
        1   -1      y   -1    -y    -y'    1     y' ;
        1   -1     -1   -y     y     1     y'   -y' ;
        1    1      y   -y    -1     y'   -y'   -1  ;
    ];

end
