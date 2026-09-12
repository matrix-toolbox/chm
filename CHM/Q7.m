function H = Q7
% ------------------------------------------------------------------------------
% 2017-03-08 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Isolated matrix Q7 taken from the doctoral thesis of F. Szollosi;
% Example 3.3.25 in https://arxiv.org/pdf/1110.5590.pdf
%
% Roots of the polynomial h(u) are given in approximate form so the numerical
% precision of H = Q7_0 is limited to 7 decimal digits. Thus the only possible
% way to get its defect is to use:
%
% >> H = Q7;
% ------------------------------------------------------------------------------
% 2022-04-16
% Q2 = sinkhorn(Q, 1e+5) provides much better signatures of CHM!
% ------------------------------------------------------------------------------

    N = 7;
    c = 1993741;
    alpha = (40169 / 3) + (50 * sqrt(c) / 3) * cos(acos(2731019453 * sqrt(c) / (c)^2) / 3 - 4 * pi / 3);

    polynomial_h = [ ...
        2054570000, ...
        4109140000, ...
        16 * alpha^2 + 9768064 * alpha - 5227993936, ...
        1956 * alpha^2 + 64132324 * alpha - 12170223176, ...
        11339 * alpha^2 + 427075897 * alpha + 1676016222, ...
        4644 * alpha^2 + 2280446676 * alpha + 8524444776, ...
        - 17074 * alpha^2 + 3269963754 * alpha + 2727593304
    ];

    rk = sort(roots(polynomial_h)); % real roots of h(u)

    z = zeros(1, 6);
    k = 1 : 6;
    z(k) = rk(k) / 2 + i * sqrt(1 - rk(k).^2 / 4)
    x = [ z(1) z(3) z(2)' z(4) z(5)' z(6)' ]; % circulant generator of Q7 core
    Q7_core = [];
    for k = 0 : N - 2
        Q7_core = [ Q7_core; circshift(x, k) ];
    end

    H = [ ones(1, N);
          ones(N - 1, 1), Q7_core ];
end
