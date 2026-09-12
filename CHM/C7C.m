function H = C7C
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated CHM of order N = 7;
% [1] J. Symb. Comput. 12, 329--336, (1991)
%
% >> H = C7C;
% ------------------------------------------------------------------------------

    a = exp(4.312838978724i); % data taken from [1]
    b = exp(1.356227956787i);
    c = exp(1.900668281165i);
    H = [
        1  1      1          1              1              1          1     ;
        1  1/a^2  a'*b'      a'*c'          a'             c/a        b/a   ;
        1  a'*b'  1/a^2/b^2  a'*c'/b^2      a'*b'*c'       c/a/b      c/a   ;
        1  a'*c'  a'*c'/b^2  1/a^2/b^2/c^2  a'/b^2/c^2     a'*b'*c'   a'    ;
        1  a'     a'*b'*c'   a'/b^2/c^2     1/a^2/b^2/c^2  1/b^2/c/a  a'*c' ;
        1  c/a    c/a/b      (a*b*c)'       1/b^2/c/a      1/a^2/b^2  a'*b' ;
        1  b/a    c/a        a'             a'*c'          a'*b'      1/a^2 ;
    ];

end
