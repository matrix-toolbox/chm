function H = S6
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Isolated BH(6, 3) found independently by G. Moorehouse and T. Tao;
%
% >> H = S6;
% ------------------------------------------------------------------------------
    H = exp(2j*pi*[
        0 0 0 0 0 0;
        0 0 1 1 2 2;
        0 1 0 2 2 1;
        0 1 2 0 1 2;
        0 2 2 1 0 1;
        0 2 1 2 1 0;
    ]/3);

end
