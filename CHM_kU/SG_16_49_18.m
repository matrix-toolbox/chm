function Y = SG_16_49_18(p)
% ------------------------------------------------------------------------------
% 2022-12-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A two-parametric self G-dual family of symmetric CHM.
% It stems from BH(16, 4) (for p1 = p2 = 0).
% ------------------------------------------------------------------------------

    Y = kron(F4(p(1)), F4(p(2)));

end

