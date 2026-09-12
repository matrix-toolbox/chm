function Y = xH_7_0_97
% ------------------------------------------------------------------------------
% 2026-09-10 Claude Opus 5
%
% isolated complex Hadamard matrix of order 7
%
% solved from xH_7_0_97_pattern_1.m
% and independently confirmed with xH_7_0_97_pattern_2.m
%
% both six-free-phase ansatzes admit a trivial attractor under random search
%    a = b = exp(-2j * pi / 3)
%    c = exp(2j * pi / 3)
%    d = 1
%    e = exp(1j * pi / 3)
%    f = -1
%
% which is exactly P_7(0) = Butson BH(7, 6) with defect 3 and #L = 6 -- NOT the intended target!
%
% the genuine L = 97 target found is monomially EQUIVALENT to Q_7 (F. Szollosi's example 3.3.25 from arXiv:1110.5590)
% it is reproduced here directly from this circulant-quotient pattern to a tighter tolerance
% no further reduction of a -- f to a closed radical form was found
% ------------------------------------------------------------------------------

    a =  0.867807015744435 + 0.496901381991173j;
    b = -0.946490532354450 - 0.322731579123256j;
    c = -0.793804522084489 - 0.608172985850420j;
    d =  0.452280350480163 - 0.891875823514429j;
    e = -0.922340249064276 + 0.386378654891868j;
    f =  0.342547937278655 + 0.939500351605117j;

    Y = [
        1    1    1    1    1    1    1;
        1    a    b    c    d    e    f;
        1    c/a  d/a  b/a  e/a  1/a  f/a;
        1    d/b  1/b  e/b  a/b  c/b  f/b;
        1    b/c  e/c  d/c  1/c  a/c  f/c;
        1    e/d  a/d  1/d  c/d  b/d  f/d;
        1    1/e  c/e  a/e  b/e  d/e  f/e;
    ];

end
