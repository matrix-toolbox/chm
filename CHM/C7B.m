function H = C7B
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated CHM of order N = 7;
% [1] J. Symb. Comput. 12, 329--336, (1991)
%
% >> H = C7B;
% ------------------------------------------------------------------------------

    d = (1j * sqrt(7) - 3) / 4;

    H = [
        1  1    1    1    1   1   1 ;
        1  d^2  d^2  d    d   1   d ;
        1  d^2  d    d^2  d   d   1 ;
        1  d    d^2  d^2  1   d   d ;
        1  d    d    1    d'  1   d';
        1  1    d    d    1   d'  d';
        1  d    1    d    d'  d'  1 ;
    ];

end
