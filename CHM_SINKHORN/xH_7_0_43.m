function Y = xH_7_0_43
% ------------------------------------------------------------------------------
% 2026-09-10 Claude Opus 5
%
% solved from xH_7_0_43_pattern_2.m
% and independently confirmed with xH_7_0_43_pattern_1.m using a different random seed
%
% isolated complex Hadamard matrix of order 7
%
% the two-free-parameter ansatz (a, c on the unit circle with b forced by the quadratic) admits TWO attractors under random search:
% - a trivial one at a = exp(2j * pi * 6 / 7) and c = exp(2j * pi / 7)
%   this is Butson BH(7, 7) = F_7 (defect 0, L = 7) -- NOT the intended target!
%
% - the genuine L = 43 target found here
%   it is monomially EQUIVALENT to C_{7C}^{(0)} [matrix-toolbox.github.io/chm/catalogue/0703.html]
%   [Haagerup, J. Symb. Comput. 12, 329-336 (1991)]
%   C7C.m --> lower-precision version
%
% no further reduction of a and c to a simple closed form (roots of unity, low-degree radicals) was found
% b follows exactly from a and c via the quadratic forced by unitarity
% ------------------------------------------------------------------------------

    a =  0.9975734668422712 - 0.06962167947049247j;
    c = -0.993359903237236  - 0.115048262222902j;

    b = ( -1 + sqrt(1 - 4*(1 + 1/a + 1/c)*(1 + a + c)) ) / (2*(1 + 1/a + 1/c));

    Y = [
        1   1             1     1        1               1            1     ;
        1   a             b     c        b*b/a           b*b          b*b/c ;
        1   a*c/b         b/c   b        c*b             b*b*b/a/c    b*b   ;
        1   a*c*c/b/b     c/b   c*c      c               c*b          b*b/a ;
        1   a*c*c/b       b/a   c*c/b    c*c             b            c     ;
        1   a/b           1/b   b/a      c/b             b/c          b     ;
        1   a*a*c*c/b/b   a/b   a*c*c/b  a*c*c/b/b       a*c/b        a     ;
    ];

end
