function H = N9
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Isolated CHM of order N = 9 found by K. Beauchamp and R. Nicoara;
% [1] https://arxiv.org/abs/math/0609076
%
% >> H = N9;
% ------------------------------------------------------------------------------
    y = (1 - 1j*sqrt(15)) / 4;

    H = [
        1,    1,    1,    1,      1,       1,       1,    1,      1;
        1,   -1,   -1,    y, y^(-1), -y^(-3),  y^(-3),    y, y^(-1);
        1,    y,  y^3, -y^2,    y^3,       y,      -1,   -y,      y;
        1,    y,    y,   -1, y^(-1),  y^(-3), -y^(-3),   -1, y^(-1);
        1,   -y, -y^2,  y^3,    y^3,      -1,       y,    y,      y;
        1,    y,  y^3,  y^3,      y,  y^(-1),  y^(-1),    y,    y^2;
        1,  y^3,    y,    y,    y^2,  y^(-1),  y^(-1),  y^3,      y;
        1, -y^2,   -y,    y,      y,      -1,       y,  y^3,    y^3;
        1,  y^3,    y,   -y,      y,       y,      -1, -y^2,    y^3;
    ];

end
