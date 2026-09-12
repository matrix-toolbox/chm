function H = C7D
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated CHM of order N = 7;
% [1] J. Symb. Comput. 12, 329--336, (1991)
%
% >> H = C7A;
% ------------------------------------------------------------------------------

    a = exp(4.312838978724i); % data taken from [1]
    b = exp(1.356227956787i);
    c = exp(1.900668281165i);
    H = [
        1  1    1        1            1            1        1  ;
        1  a^2  a*b      a*c          a            a/c      a/b;
        1  a*b  a^2*b^2  a*b^2*c      a*b*c        a*b/c    a/c;
        1  a*c  a*b^2*c  a^2*b^2*c^2  a*b^2*c^2    a*b*c    a  ;
        1  a    a*b*c    a*b^2*c^2    a^2*b^2*c^2  a*b^2*c  a*c;
        1  a/c  a*b/c    a*b*c        a*b^2*c      a^2*b^2  a*b;
        1  a/b  a/c      a            a*c          a*b      a^2;
    ];

end
