function AME_CHM36
% ------------------------------------------------------------------------------
% 2023-05-06
% ------------------------------------------------------------------------------
% take 1:
% Real HM of order 36 rotated by two permutation matrices and two unitary
% diagonal matrices with 2^nd roots of unity, does not provide (numerically)
% a 2-unitary matrix; "best" triplet of entropies is approximately
% S = 1, 0.97, 0.97...
% Maybe a matrix from another equivalence class - there are at least
% 4745357 of them - would get a better result.
% ------------------------------------------------------------------------------
% take 2
% Replace discrete diagonal matrices with continuous (complex) objects, so
% now they have q^th roots of unity with q --> infinity.
% This allows to perform continuous random walk procedure (for a moment without
% additional permutations).
% ------------------------------------------------------------------------------

    N = 36;

    H = H36;
    S_OPTIMAL = +Inf;
    MM1_OPTIMAL = zeros(N);
    MM2_OPTIMAL = zeros(N);

    while 1
        MM1 = PM(N)*DMq(N,13);
        MM2 = PM(N)*DMq(N,13);
        Y = MM1*H*MM2;
        S = abs(SL(Y) - 1) + abs(SL(R(Y)) - 1) + abs(SL(T(Y)) - 1);
        if S < S_OPTIMAL
            S_OPTIMAL = S
            SL3(MM1*H*MM2,"minimal")
            MM1_OPTIMAL = MM1;
            MM2_OPTIMAL = MM2;
        end
    end

    SL3(MM1_OPTIMAL * H36 * MM2_OPTIMAL, "compact");


end
