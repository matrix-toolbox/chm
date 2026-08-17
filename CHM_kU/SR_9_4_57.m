function Y = SR_9_4_57(p)
% ------------------------------------------------------------------------------
% 2023-01-27 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 4-parametric family of a self R-dual matrices stemming from F9_4.
% ------------------------------------------------------------------------------

    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    p4 = p(4);

    Y = F9(p); % four-parametric Fourier CHM

    r1 = rand;
    w = exp(1j*pi/9);
    DL = diag([1, 1, 1, 1, 1, 1, exp(2j*pi*r1), exp(2j*pi*r1), exp(2j*pi*r1)]);
    DR = diag([1, 1, 1, 1, w^2*exp(2j*pi*p1), w^4*exp(2j*pi*p2), 1, w^4*exp(2j*pi*p3), w^8*exp(2j*pi*p4)]);
    Y = DL * Y * DR;

end

