function Y = HH_6_4_4A
% ------------------------------------------------------------------------------
% 2023-01-15 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian BH(6, 4) of size N = 6 with d = 4 and #L = 4.
% This is a solution of "HH_6_4_34" with a = -1. It might be ~ D6.
% ------------------------------------------------------------------------------

    Y = exp(2j*pi*[
        0  0  0  0  0  0;
        0  0  2  3  2  1;
        0  2  0  2  3  1;
        0  1  2  2  0  3;
        0  2  1  0  2  3;
        0  3  3  1  1  2;
    ]/4);

end
