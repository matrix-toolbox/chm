% 20061205
% W. Bruzda, name[at]alumni.uj.edu.pl | name = w.bruzda
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/

% >> H = N11B_0
% >> abs(H .* H'), norm(H * H' - 11 * eye(11), 'fro')

function H = N11B_0

    x = (3 - i * sqrt(7))/4;
    a = -x;
    b = +x;
    y = -x^2;

    H = [
        1,  1,  1,  1,  1,  1,    1,    1,    1,    1,    1;
        1,  y,  y,  a,  b,  b,   -1,    x,   -1,    x,   -1;
        1,  y,  b,  b,  a,  y,    x,    x,   -1,   -1,   -1;
        1,  a,  b,  y,  y,  b,   -1,   -1,    x,   -1,    x;
        1,  b,  a,  y,  b,  y,    x,   -1,   -1,   -1,    x;
        1,  b,  y,  b,  y,  a,   -1,   -1,    x,    x,   -1;
        1, -1,  x, -1,  x, -1, -x/a, -x/b, -x/y, -x/y, -x/b;
        1,  x,  x, -1, -1, -1, -x/b, -x/y, -x/b, -x/y, -x/a;
        1, -1, -1,  x, -1,  x, -x/y, -x/b, -x/b, -x/a, -x/y;
        1,  x, -1, -1, -1,  x, -x/y, -x/y, -x/a, -x/b, -x/b;
        1, -1, -1,  x,  x, -1, -x/b, -x/a, -x/y, -x/b, -x/y;
    ];

end

