function [H, P1, P2, D1, D2] = BH36_6_D_2U()
% ------------------------------------------------------------------------------
% 2026-08-23
% Wojciech Bruzda | name@cft.edu.pl : name = w.bruzda
% ------------------------------------------------------------------------------
% two-unitary matrix obtained from the Dita construction by special permutations
%
% takes non-two-unitary HD = BH36_6_D() and returns H = D1 * P1 * HD * P2 * D2
% which is two-unitary:
% * rows are indexed by 6 * a + b + 1
% * columns are indexed by 6 * c + d + 1
%   with a, b, c and d in Z_6
% * all three flattenings of H are unitary
% * H is still of Butson type BH(36, 6) and it coincides with the point
%
%       BH_36((pi/3)*[0 0 0 0 0 4 5 5 4 0 0 0 0 0 0 0 0 0 0])
%
%   of the affine family of Spec. Matrices 12, 20240010 (2024)
%
% P1 and P2 are permutation matrices built from the index maps
% both are bijections by the Chinese remainder theorem, since b -> (mod(b, 2), mod(b, 3)) identifies Z_6 with Z_2 x Z_3
%
% row 6*a+b+1  ->  3*J + al + 1
%                  J = 2 * a + mod(b, 2)
%                  al = mod (b + eps(a), 3) : eps(a) = (mod(a, 3) == 2)
% col 6*c+d+1  ->  3*K + ga + 1
%                  K = 2 * c + mod(d, 2)
%                  ga = mod(-d, 3)
%
% the row map is NOT of the form (a,b) -> (perm(a), perm(b))
% the shift eps(a) inside a block depends on the leg-A index, so P1 is not a local permutation P_A kron P_B
%
% D1 and D2 are mandatory
% writing out D1 * P1 * HD * P2 * D2 entrywise gives the exact identity in Z_6
%
% log_w H(6 * a + b + 1, 6 * c + d + 1) = 4 * s(a) * al + L(J + 1, K + 1) + (be == 1) + 4 * u(c) * be + 4 * al * be
%
% where
%
% s(a) = mod(a + 1, 3)
% u(c) = mod(c, 3)
% be = mod(d, 3)
% w = exp(1j * pi / 3)
%
% the term 4 * s(a) * al depends on the row through BOTH a and al
% so it is a genuine row rescaling and no permutation of rows can produce it
% likewise (be == 1) + 4 * u(c) * be is a genuine column rescaling
%
% neither is local, which is consistent with kron(A, F) itself failing two-unitarity
%
% the column rescaling depends only on (K, be), i.e., only on the column
% so it could alternatively be absorbed into the twelve Dita blocks by calling dita_construction with B_j = F * diag(...) instead of B_j = F
% the row rescaling cannot be absorbed that way, because dita_construction lets the blocks follow the COLUMN index only
%
% >> [H _ _ _ _] = BH36_6_D_2U; summary(H), ud(H, "S", 1e-8), SL3(H)
% >> [H, P1, P2, D1, D2] = BH36_6_D_2U()
%
% ------------------------------------------------------------------------------

    w  = exp (1j * pi / 3);
    HD = BH36_6_D();

    pr = zeros (1, 36); pc = zeros (1, 36);
    dr = zeros (36, 1); dc = zeros (36, 1);
    for a = 0:5
        ea = (mod(a, 3) == 2);              % eps(a)
        sa = mod(a + 1, 3);                 % s(a)
        for b = 0:5
            J = 2*a + mod(b, 2);
            al = mod(b + ea, 3);
            pr(6 * a + b + 1) = 3 * J + al + 1;
            dr(6 * a + b + 1) = w^(mod(4 * sa * al, 6));
        end
    end
    for c = 0:5
        uc = mod(c, 3);                     % u(c)
        for d = 0:5
            K  = 2*c + mod(d, 2);
            be = mod(d, 3);
            pc(6 * c + d + 1) = 3*K + mod(-d, 3) + 1;
            dc(6 * c + d + 1) = w^(mod((be == 1) + 4 * uc * be, 6));
        end
    end

    Id = eye (36);
    P1 = Id(pr, :);                         % (P1 * HD)(i, :) = HD(pr(i), :)
    P2 = Id(:, pc);                         % (HD * P2)(:, j) = HD(:, pc(j))
    D1 = diag (dr);
    D2 = diag (dc);

    H = D1 * P1 * HD * P2 * D2;
end
