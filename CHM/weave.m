function M = weave(Aj, Bk)
% ------------------------------------------------------------------------------
% 2022-07-04 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Weaving Hadamard matrices - a construction by Robert Craigen from 1988.
% Yes: it was in nineteen eighty eight...
% -----------------------------------------------------------------------------
% Input:
%    Aj = (A1, A2, ..., An); n CHM of order m
%    Bk = (B1, B2, ..., Bm); m CHM or order n
%
% Watch out that the values of $n$ and $m$ must coincide!
% -----------------------------------------------------------------------------
% Example of call:
%
% >> M = weave({F3, F3, F3*PM(3)*DM(3), PM(3)*F3*PM(3)}, {F4(1/2), F4(1/3), F4(1/4)});
%
% where    F3 = Fourier matrix of order 3
%          F4 = 1-parameter Fourier matrix of order 4
%       PM(3) = random permutation matrix of order 3
%       DM(3) = random diagonal unitary matrix of order 3
%
% or
%
% >> x = exp(2j*pi*rand);  % random unimodular parameter for F4
% >> w = exp(2j*pi/3);     % constant 3rd root of unity for F3
%
% >> A1 = [1  1  1; 1  w  w^2; 1  w^2  w];
% >> A2 = [1  w  w; w  1  w; w  w  1];
%
% >> B1 = [1  1  1  1; 1  1 -1 -1 ; 1 -1  1 -1; 1 -1 -1  1];
% >> B2 = [1  1  i  i; 1 -1  i -i ; i  i  1  1; i -i  1 -1];
% >> B3 = [1  1  1  1; 1  1 -1 -1 ; 1 -1  x -x; 1 -1 -x  x];
%
% >> M = weave({A1, A1, A1, A2}, {B1, B2, B3});
%
% .............................................................................
% As a curiosity:
% If the free parameter in F4; x is of the form x = exp(2j*pi * n/d)
% for some non-negative integers n and d,
% then matrix M is a Butson type matrix: BH(12, q),
% where
%
%       q = lcm(12, 6*d/gcd(7*d+6*n, 5*d+6*n, 6*d))
%
% with
%       lcm = least common multiple
%       gcd = greatest common divisor
% .............................................................................
% -----------------------------------------------------------------------------

    m = size(Bk, 2);
    n = size(Aj, 2);

    % assert dimensions:
    for j = 1: n 
        assert(abs(size(Aj{j}, 1) - m) < 1e-8);
    end
    for k = 1: m 
        assert(abs(size(Bk{k}, 1) - n) < 1e-8);
    end

    % prepare matrix M by weaving Aj and Bk together
    M = zeros(m*n, m*n);
    for kk = 1 : n
    for jj = 1 : m 
        rank1submatrix = Aj{kk}(:,jj) * Bk{jj}(kk,:);
        M(1+m*(kk-1):1+m*(kk-1)+m-1, 1+n*(jj-1):1+n*(jj-1)+n-1) = rank1submatrix;
    end; end

end
