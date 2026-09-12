function H = K6_2(p)
% ------------------------------------------------------------------------------
% 2009-08-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Family of CHM depending on 2 parameters found by B. Karlsson;
% [1] https://arxiv.org/abs/0908.2555
% It looks like K6_2 allows a full domain for parameters: (p(1), p(2)) in [0, 1) x [0, 1).
%
% >> H = K6_2(rand(1, 2));
% ------------------------------------------------------------------------------

    z1 = exp(2j*pi * p(1));
    z2 = exp(2j*pi * p(2));

    Z11 = 1 - (1 - z1) * (1 - z2) / 2;
    Z12 = z2 * (1 - (1 - z1) * (1 - 1 / z2) / 2);
    Z21 = z1 * z2 * conj(Z12);
    Z22 = - z1 * z2 * conj(Z11);

    a = zeros(2);
    a(1, 1) = - Z11 * (1 + 1i * sqrt(4 / Z11 / Z11' - 1)) / 2;
    a(1, 2) = - Z12 * (1 + 1i * sqrt(4 / Z12 / Z12' - 1)) / 2;
    a(2, 1) = - Z21 * (1 - 1i * sqrt(4 / Z21 / Z21' - 1)) / 2;
    a(2, 2) = - Z22 * (1 - 1i * sqrt(4 / Z22 / Z22' - 1)) / 2;
    b = - a - [Z11 Z12; Z21 Z22];

    H = [
        1,   1,       1,       1,       1,       1;
        1,  -1,      z1,     -z1,      z1,     -z1;
        1,  z2, a(1, 1), a(1, 2), b(1, 1), b(1, 2);
        1, -z2, a(2, 1), a(2, 2), b(2, 1), b(2, 2);
        1,  z2, b(1, 1), b(1, 2), a(1, 1), a(1, 2);
        1, -z2, b(2, 1), b(2, 2), a(2, 1), a(2, 2);
    ];

end
