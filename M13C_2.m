function H = M13C_2(parameter)
% 20170330
% W. Bruzda, name[at]uj.edu.pl : name = w.bruzda
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/

% Extension of the A. LaClair's matrix M13C to the 2-parametric maximal affine family.

% >> version % 9.1.0.411655 (2016b)
% >> H = M13C_2(rand(1, 2));
% >> defect(H) % 2
% >> norm(H * H' - 13 * eye(13), 'fro')

    try
        a = parameter(1, 1);
        b = parameter(1, 2);
    catch
        warning('No valid parameter provided! Continue with default zero value.');
        a = 0;
        b = 0;
    end

    M13_C = [
        0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 5 5 5 4 3 3 2 2 2 2 0 0;
        0 2 4 2 0 1 4 1 2 4 4 0 4;
        0 2 0 0 3 2 1 3 5 5 3 5 3;
        0 2 1 4 1 4 4 5 3 1 1 5 3;
        0 2 4 1 3 4 4 1 5 1 5 3 1;
        0 0 4 2 3 5 0 4 2 0 2 2 4;
        0 5 2 2 3 2 0 5 3 3 5 5 1;
        0 4 4 0 0 1 2 4 4 2 0 2 3;
        0 0 2 4 0 3 0 2 0 2 4 3 4;
        0 4 0 4 2 5 2 0 1 4 4 2 2;
        0 4 2 2 0 5 2 2 4 5 2 4 0;
        0 2 2 4 4 1 4 4 0 4 1 2 0;
    ];

    H = exp(2 * pi * i * M13_C / 6);

    R_M13_C = [
        0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 b a 0 0 0 b 0 0 a 0 0;
        0 0 b a a a 0 b a a a a a;
        0 0 b a 0 0 0 b 0 0 a 0 0;
        0 0 b a b b 0 b b b a b b;
        0 0 b a a a 0 b a a a a a;
        0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 b a 0 0 0 b 0 0 a 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 b a b b 0 b b b a b b;
    ];

    H = H .* exp(2 * pi * i * R_M13_C);

end

