function H = HB_p_d_L(p)
% ------------------------------------------------------------------------------
% 2026-09-11 Claude Opus 5
%
% HERMITIAN Butson matrix BH(p^2, p) for every integer p >= 2, from the SYMPLECTIC form on Z_p x Z_p.
%
%     H((a,b),(c,d)) = w^(a*d - b*c),   w = exp(2j*pi/p),   a,b,c,d in BZ_p
%
% with the p^2 indices ordered lexicographically, n = p*a + b + 1.
%
% >> H = HB_p_d_L(3); % --> 9 x 9
%
% HERMITIAN
% * the symplectic form is antisymmetric, so swapping (a,b) <-> (c,d) negates the exponent:
%   H((c,d), (a,b)) = w^(c*b - d*a) = w^-(a*d - b*c) = conj(...)
% * the diagonal is a*b - b*a = 0, so H(n,n) = 1 and trace H = p^2
%   this is consistent with H^2 = p^2*I forcing eigenvalues +/- p, with multiplicities (p^2+p)/2 and (p^2-p)/2
%
% CHM
% * for two distinct indices the exponent difference is a NON-ZERO linear functional on Z_p x Z_p:
%   (a-c)*y - (b-d)*x
% * its image is the subgroup g*Z_p with g = gcd(a-c, b-d, p) < p, and every value of that image is attained the same number of times
% * summing w over a full subgroup (of size p/g > 1) gives 0, so all rows are orthogonal
%   note: this argument uses only g < p and p need NOT be prime, though BH(p^2,p) for prime p is the case usually quoted
% * it is already dephased (normalized)
%
% it is NOT symmetric (for p >= 3)
% this is the Hermitian representative of a class  whose symmetric representative is the Fourier matrix of Z_p x Z_p,
% i.e., kron(F_p, F_p) -- the two differ by the column relabelling (c,d) -> (d,-c), a pure permutation,
% so they are monomially EQUIVALENT
%
% at p = 2 the matrix is real => symmetric
%
% INVARIANTS (measured with ./defect and ./lambda, p = 2..7)
%
%     p   N=p^2    d    (p-1)(p^2-1)   #Lambda
%     2      4      3        3            2
%     3      9     16       16            3
%     4     16     57       45 (no)       4
%     5     25     96       96            5
%     6     36    259      175 (no)       6
%     7     49    288      288            7
%
%   so  #Lambda = p  for every p, while  d = (p-1)*(p^2-1)  holds for PRIME
%   p only -- at composite p the matrix is still a Hermitian BH(p^2,p) but
%   the defect is larger and no closed form is claimed here.
% ------------------------------------------------------------------------------

    if nargin < 1
        p = 3;
    endif
    if p != fix(p) || p < 2
        error('BH_pp_d_L: p must be an integer >= 2');
    endif

    N = p * p;
    w = exp(2j * pi / p);

    % indices n = p*a + b + 1  <->  (a, b)
    a = floor((0:N-1) / p);      % row/col -> first component
    b = mod((0:N-1), p);         % row/col -> second component

    % exponent matrix E(n,m) = a(n)*b(m) - b(n)*a(m)   (mod p)
    E = mod(a' * b - b' * a, p);

    H = w .^ E;

end
