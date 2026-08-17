function Y = Y_10_0_9A
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
    d = b^3;
    e = c^2;

    Y = [
        1  1  1  1  1  1  1  1  1  1  ;
        1  c  e  e  b  c  d  1  c  c  ;
        1  e  e  c  c  1  d  c  b  c  ;
        1  e  d  e  c  c  b  c  c  1  ;
        1  c  c  b  1  1  b' c' c' 1  ;
        1  1  c  c  1  b' b  1  c' c' ;
        1  c  1  c  1  c' b  c' 1  b' ;
        1  d  b  d  b' b  c  b  b' b  ;
        1  c  c  1  c' 1  b  b' 1  c' ;
        1  b  c  c  c' c' b' 1  1  1  ;
    ];

end

