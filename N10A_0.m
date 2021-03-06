function H = N10A_0
% 20061205
% W. Bruzda, name[at]uj.edu.pl : name = w.bruzda
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/

% >> version % 9.1.0.441655 (R2016b)
% >> H = N10_0
% >> abs(H .* H'), norm(H * H' - 10 * eye(10), 'fro')

    a = (i * sqrt(15) - 1) / 4;

    H = [
        1,   1,      1,      1,      1,   1,   1,   1,      1,      1;
        1,   1, a^(-2), a^(-1), a^(-2), a^2,   1, a^2,      1,      a;
        1, a^2,      1, a^(-2), a^(-2),   1,   a, a^2,      1, a^(-1);
        1,   1,      1,      1, a^(-2), a^2, a^2,   a, a^(-2), a^(-1);
        1, a^2, a^(-2),      1,      1,   a,   1, a^2, a^(-2), a^(-1);
        1,   1, a^(-1), a^(-2),      1,   1, a^2, a^2, a^(-2),      a;
        1, a^2,      1, a^(-2), a^(-1), a^2,   1,   1, a^(-2),      a;
        1, a^2, a^(-2),      1, a^(-2),   1, a^2,   1, a^(-1),      a;
        1,   a, a^(-2), a^(-2),      1, a^2, a^2,   1,      1, a^(-1);
        1,   a, a^(-1), a^(-1), a^(-1),   a,   a,   a, a^(-1),      1;
    ];

end

