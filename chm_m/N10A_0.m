% 20061205
% W. Bruzda, name[at]alumni.uj.edu.pl | name = w.bruzda
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/

% >> H = N10_0
% >> abs(H .* H'), norm(H * H' - 10 * eye(10), 'fro')

function H = N10A_0

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

