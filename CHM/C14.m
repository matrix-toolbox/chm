function H = C14
% ------------------------------------------------------------------------------
% 2023-01-08 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Isolated solution from V_N family for N = 14.
% [1] https://arxiv.org/pdf/2204.11727.pdf
%
% Y ~ C14X... : d(Y) = 0, #L(Y) = 9, it is not BH(14, q) for any q.
%
% >> H = C14;
% ------------------------------------------------------------------------------

    a = -3/4 + 1j*sqrt(7)/4;
    b = a^2;
    c = a^3;

    H = [
        1    1    1    1    1    1    1    1    1    1    1    1    1    1   ;
        1    1    1    a'   b    a    b    a    b    c    a'   a    c    a   ;
        1    1    a'   1    a    a'   a'   a    1    b    b'   b'   a    b   ;
        1    a    a'   a    b    a'   a    b    b    c    1    1    a    c   ;
        1    a    a    a    1    a    1    1    b    b    a    a    b    b   ;
        1    a'   a    1    a    1    b    b    b    a    a    a'   c    c   ;
        1    a'   c'   a'   b'   b'   1    a'   a    1    b'   c'   a    a'  ;
        1    a'   1    b'   a'   b'   a    a    1    a    1    a'   b    b   ;
        1    a'   b'   a'   a'   c'   1    b'   a    a'   c'   b'   a    1   ;
        1    b'   b'   a'   a    1    a    a'   1    b    a'   1    b    a   ;
        1    b'   c'   c'   a'   a'   b'   1    a    a'   b'   a'   1    a   ;
        1    b'   a'   c'   1    b'   a'   b'   a    a    a'   c'   a'   1   ;
        1    c'   b'   b'   b'   a'   a'   1    a    1    c'   a'   a'   a   ;
        1    c'   a'   b'   1    c'   b'   a'   a    a    a'   b'   1    a'  ;
    ];

end
