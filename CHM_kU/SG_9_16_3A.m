function Y = SG_9_16_3A(p1)
% ------------------------------------------------------------------------------
% 2022-12-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A one-parametric self G-dual family of CHM.
% It stems from BH(9, 3).
% ------------------------------------------------------------------------------

    DL = kron(diag([1, 1, exp(2j*pi*p1)]), eye(3));
    Y = DL*kron(F3, F3);

end

