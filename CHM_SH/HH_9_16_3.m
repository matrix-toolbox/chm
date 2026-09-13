function Y = HH_9_16_3()
% ------------------------------------------------------------------------------
% 2026-09-11 Claude Opus 5
%
% A Hermitian complex Hadamard matrix of size 9.
% Dimension N = 9 is a perfect square, so the "trace obstruction" does not rule it out (unlike N = 3, 5, 7, 11, ...).
%
%
% An explicit construction via the symplectic form over Z_3 x Z_3:
% index rows/columns by pairs (a, b) in Z_3^2 and set H_{(a,b), (c,d)} = omega^(ad − bc) : omega = e^(2 pi i / 3).
%
% * Hermitian: H_{(c,d), (a,b)} = omega^(cb − da) = omega^(−(ad−bc)) = conj(H_{(a,b), (c,d)})
% * unimodular: entries are cube roots of unity
% * Hadamard: for p != q the exponent difference (a−c)y − (b−d)x is a nonzero linear functional on Z_3^2
%   so it hits each of 0, 1, 2 exactly three times, moreover  3(1 + omega + omega^2) = 0
% * diagonal is all 1 (since ab−ba=0), so trace = 9, consistent with |eigenvalues| = 3 with multiplicities: n_+ = 6 and n_- = 3.
%
% Antisymmetry of the symplectic form gives H_{qp} = omega^(−(ad−bc)) = conj(H_{pq}) for free,
% and for p != q the exponent difference (a−c)y − (b−d)x is a nonzero linear functional on Z_3^2,
% so it hits 0, 1, 2 exactly three times each and 3(1+ omega + omega^2) = 0.
%
% It is permutation-equivalent to F_3 (x) F_3 with (2 5 1 4 7 0 3 6),
% exactly the relabeling (c,d) --> (d, −c) that turns the Fourier pairing ac+bd into the symplectic one ad−bc.
% So this class has both a symmetric representative (F_3 (x) F_3 itself) and a Hermitian one, which is a nice concret.
%
% It generalizes: the same formula gives a Hermitian BH(p^2, p) for every prime, so N = 4, 9, 25, 49, ...
% At p = 2 it degenerates to a real symmetric Hadamard matrix of order 4.
%
% Note: Sinkhorn under the Hermitian constraint evidently has a rough time finding this solution even though it provably exists.
% ------------------------------------------------------------------------------

    L = [
        0 0 0 0 0 0 0 0 0;
        0 0 0 2 2 2 1 1 1;
        0 0 0 1 1 1 2 2 2;
        0 1 2 0 1 2 0 1 2;
        0 1 2 2 0 1 1 2 0;
        0 1 2 1 2 0 2 0 1;
        0 2 1 0 2 1 0 2 1;
        0 2 1 2 1 0 1 0 2;
        0 2 1 1 0 2 2 1 0;
    ];

    Y = exp(2j * pi * L / 3);

end
