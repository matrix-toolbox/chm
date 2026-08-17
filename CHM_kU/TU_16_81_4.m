function Y = TU_16_81_4(p)
% ------------------------------------------------------------------------------
% 2023-01-27 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2-unitary BH(16, 4) being a result of B[8]16 x P16.
% ------------------------------------------------------------------------------

    Y = B16q4(8)*P16;
    w = exp(2j*pi/3);
    DL = diag([1 1 1 1 1 1 exp(2j*pi*p(1)) exp(2j*pi*p(1)) w w -1 -1 -1 1 -1 1]);
    DR = diag([1 1 1 1 1 i 1 i exp(2j*pi*p(2)) exp(2j*pi*p(2)) exp(2j*pi*p(2)) exp(2j*pi*p(2)) 1 i 1 i]);
    Y = DL * Y * DR;

end

