function H = P13A(p)
% ------------------------------------------------------------------------------
% 2017-07-10 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Petrescu matrix extended by F. Szollosi to 4-parametric family stemming from BH(13, 30);
% [1] https://arxiv.org/pdf/1110.5590.pdf
%
% >> H = P13A(rand(1, 4));
% ------------------------------------------------------------------------------

    a = exp(2j*pi * p(1));
    b = exp(2j*pi * p(2));
    c = exp(2j*pi * p(3));
    d = exp(2j*pi * p(4));

    w = exp(2j * pi / 5);
    o = exp(2j * pi / 3); % = omega != w :)

    p_e = acos(-cos(2*pi*p(4))/2) - 2*pi/3;
    e = exp(1j * p_e);

    assert(real(e*o) + real(d) / 2 < 1e-10)

    H = [
        o        a          b        c            -1       -a         -b       -c           w   w^2 w^3 w^4 1;
        d^2*o/a  o         -b*d*e/a -c*d*o/a/e    -d^2*o/a -1          b*d*e/a  c*d*o/a/e   w^2 w^4 w   w^3 1;
        o/b     -a*o/b/d/e -1        c*o/b/d/e    -o/b      a*o/b/d/e  o       -c*o/b/d/e   w^3 w   w^4 w^2 1;
        o/c     -a*e/c/d    b*e/c/d -1            -o/c      a*e/c/d   -b*e/c/d  o           w^4 w^3 w^2 w   1;

       -1       -a         -b       -c             o        a          b        c           w   w^2 w^3 w^4 1;
       -d^2*o/a -1          b*d*e/a  c*d*o/a/e     d^2*o/a  o         -b*d*e/a -c*d*o/a/e   w^2 w^4 w   w^3 1;
       -o/b      a*o/b/d/e  o       -c*o/b/d/e     o/b     -a*o/b/d/e -1        c*o/b/d/e   w^3 w   w^4 w^2 1;
       -o/c      a*e/c/d   -b*e/c/d  o             o/c     -a*e/c/d    b*e/c/d -1           w^4 w^3 w^2 w   1;

        w^4 w^3 w^2 w       w^4 w^3 w^2 w         1   o^2 o^2 o^2 o^2;
        w^3 w   w^4 w^2     w^3 w   w^4 w^2       o^2 1   o^2 o^2 o^2;
        w^2 w^4 w   w^3     w^2 w^4 w   w^3       o^2 o^2 1   o^2 o^2;
        w   w^2 w^3 w^4     w   w^2 w^3 w^4       o^2 o^2 o^2 1   o^2;
        1   1   1   1       1   1   1   1         o^2 o^2 o^2 o^2 1  ;
    ];

end
