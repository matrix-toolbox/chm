function Y = SR_9_16_3B
% ------------------------------------------------------------------------------
% 2023-01-27 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Example of a self R-dual matrix.
%
% It emerged as:
%
%    DL*K9_2(3)*DR
%
%    DL = diag([1, 1, 1, 1, exp(2j*pi/3), exp(4j*pi/3), 1, exp(4j*pi/3), exp(2j*pi/3)]);
%    DR = diag([1, 1, 1, 1, exp(2j*pi/3), exp(4j*pi/3), 1, exp(4j*pi/3), exp(2j*pi/3)]);
%
% and turned out to be BH(9, 3):
% ------------------------------------------------------------------------------

    Y = exp(2j*pi*[
        0   0   0   0   0   0   0   0   0;
        0   0   0   2   2   2   1   1   1;
        0   0   0   1   1   1   2   2   2;
        0   2   1   0   2   1   0   2   1;
        0   2   1   2   1   0   1   0   2;
        0   2   1   1   0   2   2   1   0;
        0   1   2   0   1   2   0   1   2;
        0   1   2   2   0   1   1   2   0;
        0   1   2   1   2   0   2   0   1;
    ]/3);

end

