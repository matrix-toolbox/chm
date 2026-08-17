function Y = TI_9_4_57A(p)
% ------------------------------------------------------------------------------
% 2022-12-29 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 4-parametric 2-unitary family of CHM.
% ------------------------------------------------------------------------------

    Y = F9(p);

    DL = diag(exp(2j*pi*[0, 1, 1, 1, 1, 0, 0, 2, 0]/3));
    DR = diag(exp(2j*pi*[0, 1, 1, 1, 0, 1, 2, 2, 1]/3));
    Y = DL * Y * DR;

end

