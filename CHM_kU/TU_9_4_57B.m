function Y = TU_9_4_57B(p)
% ------------------------------------------------------------------------------
% 2022-01-01 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Extension of T9_2(p1, p2) to 4-parametric affine family with one extra
% parameter!
%
% "exhaustive" check for defects:
%
% p1    p2    p3    p4     d
% ---------------------------
% 0  &  0  &  0  &  0  &   4* <--- generic defect
% 0  &  1  &  1  &  1  &   6
% 0  &  0  &  1  &  2  &   8
% 1  &  2  &  2  &  1  &  10
% 0  &  1  &  1  &  4  &  12
% 0  &  0  &  0  &  3  &  16
% ------------------------------------------------------------------------------
% 2023-05-05
% It looks like another representation of "TF_9_4.m"...
% ---> Removed from the official publication!
% ------------------------------------------------------------------------------

    r1 = rand;

    a = exp(2j*pi*p(1));
    b = exp(2j*pi*p(2));
    c = exp(2j*pi*p(3));
    d = exp(2j*pi*p(4));

    w=exp(2j*pi/3);

    Y = [
        1       1       1       1       1       1       1       1       1     ;
        1       w       1/w     a*w^2   a       a*w     b/w^2   b/w     b     ;
        1       1/w     w       c       c/w     c/w^2   d*w^2   d*w     d     ;
        1       1       1       w       w       w       1/w     1/w     1/w   ;
        1       w       1/w     a       a*w     a*w^2   b       b/w^2   b/w   ;
        1       1/w     w       c/w^2   c       c/w     d*w     d       d*w^2 ;
        1       1       1       1/w     1/w     1/w     w       w       w     ;
        1       w       1/w     a*w     a*w^2   a       b/w     b       b/w^2 ;
        1       1/w     w       c/w     c/w^2   c       d       d*w^2   d*w   ;
    ];

    DL = diag([1,1,1,1,w^2,w,-1j/w,-1j,-1j*w]);
    DR = diag([1,1,1,1,w,w^2,exp(2j*pi*r1),exp(2j*pi*r1)/w,exp(2j*pi*r1)/w^2]);

    Y = DL * Y * DR;

end

