function Y = BH_6_4_16
% ------------------------------------------------------------------------------
% 2022-04-10 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of size N = 6 with d = 4 and #L = 16.
% ------------------------------------------------------------------------------
% It looks like C6, but C6 is Hermitian... So, maybe it is a part of B6.
% It is NOT BH(6, q) for any q.
% Constant value of x is one of two unimodular roots {x3, x4} of monic palindromic polynomial:
%
% 1 - 2x - 2x^3 + x^4 = 0
%
% x1 = (1 - sqrt(2)*3^(1/4) + sqrt(3)) / 2
% x2 = (1 + sqrt(2)*3^(1/4) + sqrt(3)) / 2
% x3 = (1 - 1j*sqrt(2)*3^(1/4) - sqrt(3)) / 2
% x4 = (1 + 1j*sqrt(2)*3^(1/4) - sqrt(3)) / 2
%
% relations: x1x2x3x4 = 1 and x1+x2+x3+x4 = 2
% ------------------------------------------------------------------------------

    % take one root at random
    d = (1 + 1j*sqrt(2)*3^(1/4) - sqrt(3)) / 2;
    if randi(2) == 1
        d = (1 - 1j*sqrt(2)*3^(1/4) - sqrt(3)) / 2;
    end

    Y = [
        1   1     1     1     1     1  ;
        1  -1    -d    -d^2   d^2   d  ;
        1  -d    -d^3   d^4  -d^3  -d  ;
        1  -d^2   d^4  -d^4   d^2  -1  ;
        1   d^2  -d^3   d^2   1    -d' ;
        1   d    -d    -1    -d'    d' ;
    ];

end
