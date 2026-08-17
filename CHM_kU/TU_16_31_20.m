function H = TU_16_31_20(p)
% ------------------------------------------------------------------------------
% 2022-11-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% One-parametric family which yields 2-unitary CHM for particular value of "p".
%
% TO BE CONTINUED...
% ------------------------------------------------------------------------------

    c = exp(2j*pi*p);
    a = c*1j;

    H = [
        1       1       1       1       1       1       1        1       1       1        1       1      1      1      1       1   ;
        1       1j*c    1j      c       c       1j      1j*c     1      -1      -1j*c    -1j     -c     -c     -1j    -1j*c   -1   ;
        1       1j*c   -1j     -c       c       1j     -1j*c    -1       1       1j*c    -1j     -c      c      1j    -1j*c   -1   ;
        1      -1       1      -1      -1       1      -1        1       1      -1        1      -1     -1      1     -1       1   ;
        1      -1j*c   -1j/c    1       c      -1j     -1j       1/c    -1/c     1j       1j     -c     -1      1j/c   1j*c   -1   ; 
        1       c       1       c      -c      -1      -c       -1      -1      -c       -1      -c      c      1      c       1   ;
        1       c      -1      -c       c       1      -c       -1      -1      -c        1       c     -c     -1      c       1   ;
        1      -1j*c    1j/c   -1       c      -1j      1j      -1/c     1/c    -1j       1j     -c      1     -1j/c   1j*c   -1   ;
        1       1j     -1j     -1      -c      -1j/c    1j*c     1/c    -1/c    -1j*c     1j/c    c      1      1j    -1j     -1   ;
        1       c      -1      -c      -c      -1       c        1       1       c       -1      -c     -c     -1       c      1   ;
        1      -c      -1       c       c      -1      -c        1       1      -c       -1       c      c     -1      -c      1   ;
        1      -1j     -1j      1      -c       1j/c    1j*c    -1/c     1/c    -1j*c    -1j/c    c     -1      1j      1j    -1   ;
        1      -c       1/c    -1       1      -1/c     c       -1      -1       c       -1/c     1     -1      1/c    -c      1   ;
        1      -1j*c    1j     -c      -c       1j     -1j*c     1      -1       1j*c    -1j      c      c     -1j      1j*c  -1   ;
        1       1j*c    1j      c      -c      -1j     -1j*c    -1       1       1j*c     1j      c     -c     -1j     -1j*c  -1   ;
        1      -c      -1/c     1      -1       1/c     c       -1      -1       c        1/c    -1      1     -1/c    -c      1   ;
    ];

end

