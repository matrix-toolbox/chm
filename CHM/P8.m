function H = P8(p)
% ------------------------------------------------------------------------------
% 2017-12-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 3-parametric CHM of order N = 8, stemming from BH(8, 2);
% [1] F. Szollosi
%    "On quaternary complex Hadamard matrices of small orders", AMC 5 309--315 (2011)
% [2] S. Bouguezel, M.O. Ahmed, and M.N.S. Swami
%    "A new class of reciprocal-orthogonal parametric transforms", IEEE Trans. Circuits and Syst. I 56 795–805 (2009)
%
% ~~ D8_5
%
% >> H = P8(rand(1, 3));
% ------------------------------------------------------------------------------

    a = exp(2j*pi * p(1));
    b = exp(2j*pi * p(2));
    c = exp(2j*pi * p(3));
    H = [
        1  1  1  1  1  1  1  1;
        1 -a  b -c  c -b  a -1;
        1  a -b -c  c  b -a -1;
        1 -1 -1  1  1 -1 -1  1;
        1  a  b  c -c -b -a -1;
        1 -1  1 -1 -1  1 -1  1;
        1  1 -1 -1 -1 -1  1  1;
        1 -a -b  c -c  b  a -1;
    ];

    % it may be derived from BH(8, 2)
    %
    %   0   0   0   0   0   0   0   0
    %   0   1   0   1   0   1   0   1
    %   0   0   1   1   0   0   1   1
    %   0   1   1   0   0   1   1   0
    %   0   0   0   0   1   1   1   1
    %   0   1   0   1   1   0   1   0
    %   0   0   1   1   1   1   0   0
    %   0   1   1   0   1   0   0   1

end
