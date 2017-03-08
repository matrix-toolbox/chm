% 20061205
% W. Bruzda, name[at]alumni.uj.edu.pl | name = w.bruzda
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/

% >> H = F5_0
% >> abs(H .* H'),  norm(H * H' - 5 * eye(5), 'fro')

function H = F5_0

    w = exp(2 * pi * i / 5);

    H = [
        1,   1,   1,   1,   1;
        1,   w, w^2, w^3, w^4;
        1, w^2, w^4,   w, w^3;
        1, w^3,   w, w^4, w^2;
        1, w^4, w^3, w^2, w
    ];

end

