function Y = Y_10_0_9B
% -----------------------------------------------------------------------------
% 2023-01-07 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% An isolated matrix found by the Sinkhorn algorithm.
% It is N9...
% -----------------------------------------------------------------------------

    b = -1/4 -1j*sqrt(15)/4;
    c = b^2;

    Y = [
        1  1  1  1  1  1  1  1  1  1  ;
        1  1  1  1  b  b  c  c  c  c  ;
        1  1  c  c  c  b  1  b  c  1  ;
        1  c  1  c  c  b  c  1  1  b  ;
        1  c  c  b  1  b  1  c  1  c  ;
        1  b  b  b' b' 1  b  b' b  b' ;
        1  1  c' c' 1  b' c' 1  b' c' ;
        1  c' 1  c' 1  b' b' c' c' 1  ;
        1  b' c' 1  c' b' c' c' 1  1  ;
        1  c' b' 1  c' b' 1  1  c' c' ;
    ];

end

