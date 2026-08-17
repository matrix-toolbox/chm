function Y = Y_10_0_6
% -----------------------------------------------------------------------------
% 2023-01-07 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% An isolated and symmetric BH(10, 6) in the LOG-form
% found by the Sinkhorn algorithm.
%
% [1] https://arxiv.org/abs/2204.11727
% -----------------------------------------------------------------------------

    Y = exp(2j*pi*[
        0  0  0  0  0  0  0  0  0  0;
        0  0  3  3  3  3  1  1  5  5;
        0  3  3  4  5  1  5  1  3  1;
        0  3  4  3  1  5  1  5  1  3;
        0  3  5  1  3  2  3  5  5  1;
        0  3  1  5  2  3  5  3  1  5;
        0  1  5  1  3  5  5  2  3  3;
        0  1  1  5  5  3  2  5  3  3;
        0  5  3  1  5  1  3  3  1  4;
        0  5  1  3  1  5  3  3  4  1;
    ]/6);

end

    % % alternative (original) form
    %
    % w = exp(2j*pi/6);
    %
    % Y = [
    %     1   1     1     1     1      1     1     1     1     1   ;
    %     1   1    -1    -1    -1     -1     w     w     1/w   1/w ;
    %     1  -1    -1    -w     1/w    w     1/w   w    -1     w   ;
    %     1  -1    -w    -1     w      1/w   w     1/w   w    -1   ;
    %     1  -1     1/w   w    -1     -1/w  -1     1/w   1/w   w   ;
    %     1  -1     w     1/w  -1/w   -1     1/w  -1     w     1/w ;
    %     1   w     1/w   w    -1      1/w   1/w  -1/w  -1    -1   ;
    %     1   w     w     1/w   1/w   -1    -1/w   1/w  -1    -1   ;
    %     1   1/w  -1     w     1/w    w    -1    -1     w    -w   ;
    %     1   1/w   w    -1     w      1/w  -1    -1    -w     w   ;
    % ];

