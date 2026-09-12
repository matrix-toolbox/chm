function H = Q11X(k)
% ------------------------------------------------------------------------------
% 2017-03-08 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated matrix Q11 taken from the doctoral thesis of F. Szollosi;
% Example C.2.11 in https://arxiv.org/pdf/1110.5590.pdf
%
% Symmetric isolated matrix with d = 0 and #L = 63.
% There are 2 pairs of complex conjugate roots of the polynomial h(a)
% given in approximate form... See Q7_0.m for comparison.
%
% >>
% >> k = 1 % 2 3 or 4
% >> H = Q11X_0(k);
% ------------------------------------------------------------------------------

    N = 11;
    gamma = sqrt(396 + 18 * (8468 + 12 * sqrt(12309))^(1/3) + 3 * (1829088 - 2592 * sqrt(12309))^(1/3)) / 3;
    polynomial_h = [ ...
        1008, ...
        -(336 * gamma - 1344), ...
        -(gamma^5 - 132 * gamma^3 - 56 * gamma^2 + 1088 * gamma + 1232), ...
        -(336 * gamma - 1344), ...
        1008
    ];
    s = roots(polynomial_h);
    try
        a = s(SIGMA); % polynomial has two (complex conjugate) pairs of roots
    catch
        warning("No valid SIGMA = 1 2 3 or 4 parameter provided! Continue with default value: SIGMA = 1");
        a = s(1);
    end

    b_plus_c = - (30*a^7 + 53*a^6 - 83*a^5 - 56*a^4 + 117*a^3 -176*a^2 + 103*a + 59) / 48;
    t = fsolve(@(x) bc_decomposition(x, b_plus_c), [0 0 0 0]);
    b = t(1) + i * t(2);
    c = t(3) + i * t(4);
    x = [ a b b' c c' a' c' c b' b ]; % circulant generator of Q11_SIGMA core
    Q11_core = [];
    for k = 0 : N - 2
        Q11_core = [ Q11_core; circshift(x, k) ];
    end

    H = [ ones(1, N);
          ones(N - 1, 1), Q11_core ];
end

% ------------------------------------------------------------------------------

function F = bc_decomposition(x, b_plus_c)
    F = [
        x(1) * x(1) + x(2) * x(2) - 1;
        x(3) * x(3) + x(4) * x(4) - 1;
        x(1) + x(3) - real(b_plus_c);
        x(2) + x(4) - imag(b_plus_c);
   ];
end
