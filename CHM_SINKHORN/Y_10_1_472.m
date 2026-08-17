function Y = Y_10_1_472(p1)
% -----------------------------------------------------------------------------
% 2022-08-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% A 1-parametric non-affine family of CHM of order N = 10 found by the Sinkhorn algorithm.
% It has generic defect = 1 and #L = 472.
% Phase p1 does not run over full circle!
%
% [1] https://arxiv.org/abs/2204.11727
% -----------------------------------------------------------------------------

    I = 1j;

    a = exp(2j * pi * p1);
    b = -I - (2 * a) / (1 + a * (1 + I + a));

    zeta1 = sqrt(3 + 4 * a + 6 * a^2 + 4 * a^3 + 3 * a^4);

    c = -((1 + 2 * a + a^2 + I * zeta1) / (2 * (1 + a + a^2)));

    zeta2 =  1 - 2 * a^2 - 4 * a^3 - 5 * a^4 - 4 * a^5 - 2 * a^6 - I * zeta1 - 2 * I * a * zeta1 - I * a^2 * zeta1;
    zeta3 =  a + 4 * a^2 + 4 * a^3 + 3 * a^4 + a^5 + a^6 + I * a * zeta1 - I * a^3 * zeta1 - I * a^4 * zeta1;
    zeta4 = -2 - 3 * a - 3 * a^2 + a^4 + a^5 + I * a * zeta1 + I * a^2 * zeta1 + I * a^3 * zeta1;

    d = (zeta2 - sqrt(zeta2^2 - 4 * zeta3 * zeta4)) / 2 / zeta3;
    e = -((I * (1 + a * ((1 + I) + a)) * d * (-I + a * c)) / (1 + a * ((1 + I) + a + I * (1 + a * ((1 - I) + a)) * d) + d * (-I + a * (d + a * (I + (1 + I) * d + a * ((1 + I) + d)))) * c));
    f = -(((1 + a * ((1 - I) + a)) * (1 + d * ((1 - I) + d))) / ((1 + a * ((1 + I) + a)) * (1 + d * ((1 + I) + d))));

    Y = [
        1    1        1         1            1         1        1         1             1           1;
        1   -1        I         I           -I        -I       -I*a*c     I*a*c        -I/a/c       I/a/c;
        1    I        f        -f            b         I*b      a*b       I*a           b/a         I/a;
        1    I       -f         f           -I*b      -b       -b*c       I*c          -b/c         I/c;
        1   -I        f/b      -I*f/b        I*f      -I*f      1/d      -I*f/b/d       d          -I*f*d/b;
        1   -I        I*f/b    -f/b         -I*f       I*f      a*d*c     I*a*f*d*c/b   1/a/d/c     I*f/a/b/d/c;
        1   -I*a*c    f/b/d    -a*f*d*c/b    a         c        I*a*c    -I*a*e*c       I           e;
        1    I*a*c    I/d       I*a*d*c     -I*a*b     I*b*c   -a*c/e    -a*c           I/e        -1;
        1   -I/a/c    f*d/b    -f/a/b/d/c    1/a       1/c      I        -e             I/a/c       I*e/a/c;
        1    I/a/c    I*d       I/a/d/c     -I*b/a     I*b/c   -I/e      -1             1/a/e/c    -1/a/c;
    ];

end
