function Y = SH_8_9_10(p)
% ------------------------------------------------------------------------------
% 2023-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric 1-parameter family of CHM of order N = 8
% with generic d = 9 and #L = 10.
% It stems from BH(8, 2).
% It was originally found as a solution of "SH_8_3_142.m".
% ------------------------------------------------------------------------------

    c = exp(2j*pi* p(1));

    Y = [
        1    1    1    1    1     1    1     1  ;
        1    1   -1   -c    c     c   -c    -1  ;
        1   -1    1    1   -1     1   -1    -1  ;
        1   -c    1   -1    c    -1    c    -c  ;
        1    c   -1    c   -c^2  -c    c^2  -c  ;
        1    c    1   -1   -c    -1   -c     c  ;
        1   -c   -1    c    c^2  -c   -c^2   c  ;
        1   -1   -1   -c   -c     c    c     1  ;
    ];

end

