function Y = xH_10_1_L(t)
% ------------------------------------------------------------------------------
% 2026-09-10 Claude Opus 5
%
% one-parameter affine family of complex Hadamard matrices of order 10 and defect d = 1
% recovered from numerical data
%
% >> H = xH_10_1_L(t); % a = exp(2j * pi * t) is the single free real phase
%
% all other constants: b, c, d, e, f and the auxiliary zeta's
% are forced by unitarity and follow from a via the closed-form radicals below
%
% t in (0, 0.25) or (0.369, 0.63)
%
% generic #L = 472
%
% t = 0    =>  #L = 208
% t = 0.1  =>  #L = 472
% t = 0.2  =>  #L = 396
% t = 0.24 =>  #L = 472
% t = 0.25 =>  ...
% t = 0.3  =>  NOT CHM
% t = 0.4  =>  #L = 396
% t = 0.5  =>  #L = 8
                d = 8
                class BH(10, 8)
% t = 0.55 =>  #L = 472
% t = 0.6  =>  #L = 396
% ------------------------------------------------------------------------------

    if nargin == 0
        t = 0.5;
    endif

    a = exp(2j*pi*t);

    c = -1j - 2*a / (1 + a*((1 + 1j) + a));

    zeta  = -sqrt(3 + 4*a + 6*a^2 + 4*a^3 + 3*a^4);
    zeta2 = 1 - 2*a^2 - 4*a^3 - 5*a^4 - 4*a^5 - 2*a^6 - 1j*zeta - 2j*a*zeta - 1j*a^2*zeta;
    zeta3 = a + 4*a^2 + 4*a^3 + 3*a^4 + a^5 + a^6 + 1j*a*zeta - 1j*a^3*zeta - 1j*a^4*zeta;
    zeta4 = -2 - 3*a - 3*a^2 + a^4 + a^5 + 1j*a*zeta + 1j*a^2*zeta + 1j*a^3*zeta;

    e = (zeta2 - sqrt(zeta2^2 - 4*zeta4*zeta3)) / (2*zeta3);
    b = -((1 + a*((1 - 1j) + a)) * (1 + e*((1 - 1j) + e))) / ((1 + a*((1 + 1j) + a)) * (1 + e*((1 + 1j) + e)));
    f = -(1 + 2*a + a^2 + 1j*zeta) / (2*(1 + a + a^2));

    p = 1j*(1 + a*((1 + 1j) + a)) * e * (-1j + a*f);
    q = 1 + a*((1 + 1j) + a + 1j*(1 + a*((1 - 1j) + a))*e) + e*(-1j + a*(e + a*(1j + (1 + 1j)*e + a*((1 + 1j) + e)))) * f;
    d = -p/q;

    Y = [
        1  1       1         1           1       1        1        1              1         1            ;
        1 -1       1j        1j         -1j     -1j      -1j*a*f   1j*a*f        -1j/a/f    1j/a/f       ;
        1  1j      b        -b           c       1j*c     a*c      1j*a           c/a       1j/a         ;
        1  1j     -b         b          -1j*c   -c       -c*f      1j*f          -c/f       1j/f         ;
        1 -1j      b/c      -1j*b/c      1j*b   -1j*b     1/e     -1j*(b/c)/e     e        -1j*b*e/c     ;
        1 -1j      1j*b/c   -b/c        -1j*b    1j*b     a*e*f    1j*a*b*e*f/c   1/a/e/f   1j*b/a/c/e/f ;
        1 -1j*a*f  (b/c)/e  -a*b*e*f/c   a       f        1j*a*f  -1j*a*d*f       1j        d            ;
        1  1j*a*f  1j/e      1j*a*e*f   -1j*a*c  1j*c*f  -a*f/d   -a*f            1j/d     -1            ;
        1 -1j/a/f  b*e/c    -b/a/c/e/f   1/a     1/f      1j      -d              1j/a/f    1j*d/a/f     ;
        1  1j/a/f  1j*e      1j/a/e/f   -1j*c/a  1j*c/f  -1j/d    -1              1/a/d/f  -1/a/f        ;
    ];

end
