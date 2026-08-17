function Y = SG_9_16_3B()
% ------------------------------------------------------------------------------
% 2022-12-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Self G-dual matrix derived from the Karlsson's family.
% ------------------------------------------------------------------------------

    w = exp(2j*pi/3);
    DL = diag([1 1 1 1 w^2 w w w^2 1]);

    Y = DL * K9_2z(3) * DL;

end

