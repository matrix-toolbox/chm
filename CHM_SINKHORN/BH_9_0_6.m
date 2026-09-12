function Y=BH_9_0_6
% ------------------------------------------------------------------------------
% 2023-01-06 Wojciech Bruzda, name[at]uj.edu.pl : name = w.bruzda, https://matrix-toolbox.github.io/chm
% ------------------------------------------------------------------------------
% An isolated BH(9, 6) matrix found by the Sinkhorn algorithm.
%
% [1] https://arxiv.org/abs/2204.11727
% ------------------------------------------------------------------------------

    w = exp(2j*pi/6);

    Y = [
        1   1   1   1   1   1   1   1   1  ;
        1   1   1   w  -w  -w  -w' -1  -1  ;
        1   1  -w  -1  -w' -w' -w'  w'  w' ;
        1  -w' -w  -1   1  -w   1  -1   w  ;
        1  -w' -w   1  -w'  1  -w  -w' -w  ;
        1  -w' -w'  w' -w  -w'  1   w' -1  ;
        1  -w  -w' -w  -w   1  -w' -w'  1  ;
        1  -w   1  -1   1  -w' -w   w  -1  ;
        1  -w  -w'  w  -w' -w  -w   w'  w  ;
    ];

end
