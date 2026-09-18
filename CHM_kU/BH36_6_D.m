function [H, A, F]= BH36_6_D()
% ------------------------------------------------------------------------------
% 2026-08-23
% Wojciech Bruzda | name@cft.edu.pl : name = w.bruzda
%
% DC-observation due to Tuomo Valtonen (more to come soon!)
% ------------------------------------------------------------------------------
% H = dita_construction (A, F, F, ..., F) = kron (A, F)
%
% A : the BH(12, 6) core given below, A = w.^L with w = exp(1j * pi / 3)
% F : the isolated Fourier matrix of order three
%
% H is a CHM of order 36 by Dita's theorem
%
% H IS NOT TWO-UNITARY
% splitting the row index of a tensor-product matrix as (a, b)
% with a, b in Z_6, the reshuffled unitarity condition would demand
%
% sum_{b2, d2 in Z_2} A((a, b2), (c, d2)) * conj (A((a', b2), (c', d2)))
% = 4 * delta(a, a') * delta(c, c')
%
% i.e., 36 pairwise orthogonal vectors in C^4 --> impossible
%
% two-unitarity has to be restored by hand, which done in BH_36_6_D_2U
% ------------------------------------------------------------------------------

    w = exp(1i * pi / 3);

    L = [
        5 3 2 2 2 4 5 3 4 4 4 2
        5 1 1 1 0 4 2 4 0 0 1 3
        4 1 5 4 3 4 0 3 1 2 3 4
        0 3 0 1 3 2 5 2 1 0 0 5
        1 3 0 4 4 4 3 1 0 4 2 2
        5 3 3 5 0 0 0 2 0 2 1 1
        5 0 4 1 4 5 5 0 2 5 2 1
        5 4 3 0 4 3 2 1 4 1 3 4
        0 0 1 5 3 1 4 4 5 1 3 1
        2 2 4 0 3 5 3 3 3 1 0 2
        3 4 0 1 2 5 1 0 0 1 4 1
        3 2 3 2 4 1 2 3 0 5 3 0
    ];
    A = w.^L;

    F = F3 ();

    % check the format!
    % the permutations of BH_36_D_2U assume
    % F(al + 1, be + 1) = w3^(al * be) with w3 = exp(2i * pi / 3)
    % if F3 returns the conjugate convention --> flip it here so the rest stays valid
    if abs (F(2, 2) - exp(-2j * pi / 3)) < 1e-9
        F = conj (F);
    elseif abs (F(2, 2) - exp(2j * pi / 3)) > 1e-9
        error ("F3 is not a Fourier matrix in either convention!");
    end

    % twelve identical order-three blocks --> the plain tensor product kron(A, F)
    H = dita_construction(A, F, F, F, F, F, F, F, F, F, F, F, F);
end
