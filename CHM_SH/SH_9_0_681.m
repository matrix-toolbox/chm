function Y = SH_9_0_681
% ------------------------------------------------------------------------------
% 2026-09-13 Wojciech Bruzda, name[at]uj.edu.pl : name = w.bruzda, https://matrix-toolbox.github.io/chm
% ------------------------------------------------------------------------------
% Symmetric CHM of order N = 9 with d = 0 and #L = 681, built from twenty numbers.
%
% Dephased at the right entry the 8x8 core takes only 20 distinct values out of
% 64. It is symmetric, and invariant under simultaneous row and column
% permutations forming a group of order 2, generated on the core indices 1..8 by the involution
%
%     sigma = (1 8)(2 3)(4 6)(5 7)
%
% so each orbit on the 64 core positions carries one value.
% Brute force over all 8! permutations confirms the group is exactly of
% order 2, and Burnside then gives the orbit count (64 + 8 + 0 + 8)/4 = 20,
% against the 36 a generic symmetric 8x8 core would need. All eight
% #L = 681 matrices found in a random search share this structure, so it is
% a property of the stratum rather than an accident.
%
% Unlike the N = 11, #L = 139 case the core is not circulant and no dephasing
% makes it so, which is why the reduction stops at 20 free phases.
%
% >> Y = SH_9_0_681;
% ------------------------------------------------------------------------------

    p = [
        0.91334665244472535,     % a
        0.67077265916958873,     % b
        0.36729110515828872,     % c
        0.61955220061331351,     % d
        0.22160654567837285,     % e
        0.94538716674496936,     % f
        0.27235027038642062,     % g
        0.54837092842665447,     % h
        0.87358954846271208,     % i
        0.62467885889676500,     % j
        0.20904155942955360,     % k
        0.11851912375211690,     % l
        0.33543226398371279,     % m
        0.69795810272968750,     % n
        0.47529333649141703,     % o
        0.53472276363079319,     % p
        0.04077821601779233,     % q
        0.75683996733724601,     % r
        0.81245074359153779,     % s
        0.42846041648672795;     % t
    ];

    % which of the 20 phases each entry of the 8x8 core uses
    C = [
         1  2  3  4  5  6  7  8;
         2  9 10 11 12 13 14  3;
         3 10  9 13 14 11 12  2;
         4 11 13 15 16 17 18  6;
         5 12 14 16 19 18 20  7;
         6 13 11 17 18 15 16  4;
         7 14 12 18 20 16 19  5;
         8  3  2  6  7  4  5  1;
    ];

    Y = ones(9);
    Y(2:9, 2:9) = exp(2j * pi * p(C));

end
