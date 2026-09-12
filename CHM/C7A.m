function H = C7A
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated CHM of order N = 7;
% [1] J. Symb. Comput. 12, 329--336, (1991)
%
% >> H = C7A;
% ------------------------------------------------------------------------------

    d = (1j * sqrt(7) - 3) / 4;

    H = [
        1 1     1      1      1  1  1 ;
        1 1/d^2 1/d^2  d'     d' 1  d';
        1 1/d^2 d'     1/d^2  d' d' 1 ;
        1 d'    1/d^2  1/d^2  1  d' d';
        1 d'    d'     1      d  1  d ;
        1 1     d'     d'     1  d  d ;
        1 d'    1      d'     d  d  1 ;
    ];

    % k_0 = (5 - sqrt(21)) / 4 - i * (3 / 2) * sqrt((sqrt(21) - 3) / 2); % there is a typo in the formulas in Haagerup's paper!
    % k_1 = (3 - sqrt(21)) / 4 - i * 2 * sqrt((sqrt(21) - 3) / 2);
    % k_2 = i * sqrt((sqrt(21) + 3) / 2);
    % A = k_0 + k_1 * cos(2 * pi / 7) + k_2 * sin(2 * pi / 7);
    % B = k_0 + k_1 * cos(8 * pi / 7) + k_2 * sin(8 * pi / 7);
    % C = k_0 + k_1 * cos(4 * pi / 7) + k_2 * sin(4 * pi / 7);
    % a = A
    % b = B / A
    % c = C / B

end
