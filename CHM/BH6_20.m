function H = BH6_20
% ------------------------------------------------------------------------------
% 2021-04-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Example of B(6, 20) with defect = 4 and #L = 12.
%
% Quest: Is it possible to symmetrize it?
%
% >> H = BH6_20;
% ------------------------------------------------------------------------------

    H = exp(2j*pi*[
        0   0   0   0    0   0;
        0   0   5  15   10  10;
        0   5  15  10    2  12;
        0  15  10   5    7  17;
        0  10  17   2   12   7;
        0  10   7  12   17   2;
    ]/20);

end
