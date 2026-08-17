function Y = FH_7_0_14
% ------------------------------------------------------------------------------
% 2022-04-09 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated BH(7, 14) ~ F7; found by the Sinkhorn algorithm.
% ------------------------------------------------------------------------------

    Y = exp(2j*pi*[
        7   0   0   0   0   0   0;
        0   5   9  13   1  11   3;
        0   9   5   1  13   3  11;
        0  13   1   3  11   9   5;
        0   1  13  11   3   5   9;
        0  11   3   9   5  13   1;
        0   3  11   5   9   1  13;
    ]/14);

end

