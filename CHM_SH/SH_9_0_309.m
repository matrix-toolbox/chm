function Y = SH_9_0_309
% ------------------------------------------------------------------------------
% 2026-09-13 Wojciech Bruzda, name[at]uj.edu.pl : name = w.bruzda, https://matrix-toolbox.github.io/chm
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 9 with d = 0 and #L = 309, built from thirteen numbers.
%
% Dephased at the right entry the 8x8 core takes only 13 distinct values out of
% 64. It is symmetric, and invariant under simultaneous row and column
% permutations forming a Klein four-group of order 4, generated on the core indices 0..7 by
%
%     (0 1)(2 7)(3 6)(4 5)   and   (0 4)(1 5)(3 6)
%
% so each orbit on the 64 core positions carries one value.
% Brute force over all 8! permutations confirms the group is exactly this.
%
% >> Y = SH_9_0_309;
% ------------------------------------------------------------------------------

    p = [
        0.03763910451805388,     % a
        0.73605334227920960,     % b
        0.48436077927330257,     % c
        0.36376036880564716,     % d
        0.26053207127440914,     % e
        0.32593926568088799,     % f
        0.76811216934591231,     % g
        0.77844708109786187,     % h
        0.11975423930448972,     % i
        0.21847570251379769,     % j
        0.66871559691051996,     % k
        0.83030726203713023,     % l
        0.60420156584632456;     % m
    ];

    % which of the 13 phases each entry of the 8x8 core uses
    C = [
         1  2  3  4  5  6  7  8;
         2  1  8  7  6  5  4  3;
         3  8  9 10  3  8 10 11;
         4  7 10 12  7  4 13 10;
         5  6  3  7  1  2  4  8;
         6  5  8  4  2  1  7  3;
         7  4 10 13  4  7 12 10;
         8  3 11 10  8  3 10  9;
    ];

    Y = ones(9);
    Y(2:9, 2:9) = exp(2j * pi * p(C));

end
