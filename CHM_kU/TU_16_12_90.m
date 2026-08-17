function Y = TU_16_12_90(p)
% ------------------------------------------------------------------------------
% 2022-12-30 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% When matrix, being an output from >> sfc(exp(2j*pi*A), "VERBOSE");
% is given in the "symbolic" from with variables from a to (say) q,
% assume that all of them can be used as free affine parameters.
% Usually this does not work, because
%    nh(Z) != 0
%    n1(Z) = 0
% So, try to replace one variable with something that depends on another one, eg. : b --> 1/c/c (get rid of b).
% Check again nh(Z) and n1(Z).
% Repeat this procedure until Z(a, b, ...) becomes CHM.
% ------------------------------------------------------------------------------
%
% 3-parameter 2-unitary CHM family.
%
%
% Check again this table: all parameters should have the form integer/16 or integer/32...
%
% p1    p2    p3    d    #L   q : BH(16, q)
%-------------------------------------------
% r     r     r     12   -    -
% 1/5   1/4   1/4   12   20   20
% 0     0     1/8   14   8    8
% 1/8   0     0     16   8    8
% 11/16 13/16 12/16 17   16   16
% 0     3/8   1/8   19   8    8
% 0     1/8   1/4   20   8    8
% 15/16 13/16 12/16 21   16   16
% 1/8   1/8   0     25   8    8
% 0     1/8   1/8   29   8    8
% 2/3   3/4   2/3   31   12   12
% 1/4   1/3   1/2   35   12   12
% 1/4   1/4   1/4   37   4    4
% 1/8   0     1/8   53   8    8
% 1/4   0     1/4   57   4    4
% 1/2   1/3   1/2   70   6    6
% 0     r     0     70   -    -
% 0     3/4   1/2   71   4    4
% 0     0     0     105  2    2
%------------------------------------------
%
%
%
% "exhaustive" check for Z(p1/16,p2/16,p3/16), * marks generic value of defect
% p1 &  p2 &  p3 &  d
%--------------------
% 1  &  0  &  0  &  12 *
% 0  &  0  &  1  &  14
% 1  &  0  &  6  &  16
% 1  &  2  &  7  &  17
% 0  &  1  &  5  &  19
% 0  &  1  &  4  &  20
% 1  &  2  &  3  &  21
% 0  &  1  &  1  &  25
% 0  &  2  &  2  &  29
% 1  &  4  &  1  &  31
% 4  &  1  &  0  &  35
% 0  &  0  &  4  &  37
% 1  &  0  &  1  &  53
% 4  &  0  &  4  &  57
% 0  &  1  &  0  &  70
% 0  &  4  &  0  &  71
% 0  &  0  &  0  &  105
% ------------------------------------------------------------------------------

    a = exp(2j*pi*p(1));
    b = exp(2j*pi*p(2));
    c = exp(2j*pi*p(3));

    Y = [
    1  1          1      1          1          1        1     1        1      1        1        1        1        1       1           1        ;
    1  b*a/c     -b*a/c  b*a/c      b         -b       -1    -b        b*a/c -b*a/c   -c/a     -b*a/c    c/a      b       -b          b        ;
    1 -1          c/a/b -1          c^2/b/a^2  a        1    -a       -c/a/b -c/a^2   -c/a      c/a^2   -c/a      c/a     -c^2/b/a^2  c/a      ;
    1  1          1      1         -1         -1        1    -1        1      1       -1        1       -1       -1       -1         -1        ;
    1  a^2/c     -a/c   -a^2/c      1         -c/a     -a/c  -c/a     -1      a/c      c/a      a/c     -1       -c/a^2    c/a        c/a^2    ;
    1 -a^2/c     -c/a/b  a^2/c     -c^2/b/a^2  a       -1    -a        c/a/b -c/a^2    c/a      c/a^2   -c/a      c^2/a^3  c^2/b/a^2 -c^2/a^3  ;
    1 -b*a/c      b*a/c -b*a/c      b          b*a^2/c -1    -b*a^2/c  b*a/c  b/a     -c/a     -b/a      c/a     -b        b         -b        ;
    1  a^2/c      a/c   -a^2/c     -1         -c/a      a/c  -c/a     -1     -a/c      c/a     -a/c      1        c/a^2    c/a       -c/a^2    ;
    1 -1          c/a   -1         -c^2/a^2   -a       -c/a   a        1      c^2/a^3  c/a     -c^2/a^3 -c^2/a^2  c^2/a^2 -c/a        c^2/a^2  ;
    1  1          c/a/b  1          c^2/b/a^2  c/a     -1     c/a     -c/a/b -1        c/a     -1       -c/a     -c/a     -c^2/b/a^2 -c/a      ;
    1  b*a^3/c^2 -b*a/c -b*a^3/c^2 -b          b       -1     b       -b*a/c  b*a/c   -c/a      b*a/c    c/a      b*c/a^2 -b         -b*c/a^2  ;
    1 -1         -c/a   -1          c^2/a^2   -a        c/a   a        1     -c^2/a^3  c/a      c^2/a^3  c^2/a^2 -c^2/a^2 -c/a       -c^2/a^2  ;
    1 -a^2/c     -1      a^2/c     -c^2/a^2    c        1    -c       -1      c/a^2    c^2/a^2 -c/a^2    c^2/a^2 -c^3/a^4 -c^2/a^2    c^3/a^4  ;
    1 -b*a^3/c^2  b*a/c  b*a^3/c^2 -b         -b*a^2/c -1     b*a^2/c -b*a/c -b/a     -c/a      b/a      c/a     -b*c/a^2  b          b*c/a^2  ;
    1 a^2/c      -c/a/b -a^2/c     -c^2/b/a^2  c/a      1     c/a      c/a/b -1       -c/a     -1       -c/a     -c^2/a^3  c^2/b/a^2  c^2/a^3  ;
    1 -a^2/c     -1      a^2/c      c^2/a^2    -c       1     c       -1      c/a^2   -c^2/a^2 -c/a^2   -c^2/a^2  c^3/a^4  c^2/a^2   -c^3/a^4  ;
    ];


    DL = diag([1,1,1,1,1,a*b/c, -c/a/b,1,-1,b,c^2/b/a^2,1,-1,c/a,-c/a,1]);
    DR = diag([1,1,1,1,1,-1,-c/a,-1,c/a/b,-a/b,1/b,a/b,-c/a^2,-1,c/a^2,1]);
    Y = DL * Y * DR;

end

