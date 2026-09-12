function H = F6(p)
% ------------------------------------------------------------------------------
% 2016-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Fourier CHM depending on 2 parameters;
%
% >> H = F6(rand(1, 2));
% ------------------------------------------------------------------------------

    R = [
        0 0 0 0 0 0;
        0 p(1) p(2) 0 p(1) p(2);
        0 0 0 0 0 0;
        0 p(1) p(2) 0 p(1) p(2);
        0 0 0 0 0 0;
        0 p(1) p(2) 0 p(1) p(2);
    ];

    H = FN(6) .* exp(2j * pi* R);

end
