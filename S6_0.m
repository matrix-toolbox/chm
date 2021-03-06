function H = S6_0
% 20061205
% W. Bruzda, name[at]uj.edu.pl : name = w.bruzda
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/

% >> version % 9.1.0.441655 (R2016b)
% >> H = S6_0
% >> abs(H .* H'), norm(H * H' - 6 * eye(6), 'fro')

    w = exp(2 * pi * i / 3);

    H = [
        1,       1,       1,       1,       1,     1;
        1,       1,       w,       w,     w^2,   w^2;
        1,       w,       1,     w^2,     w^2,     w;
        1,       w,     w^2,       1,       w,   w^2;
        1,     w^2,     w^2,       w,       1,     w;
        1,     w^2,       w,     w^2,       w,     1;
    ];

end

