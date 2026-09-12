function H = G10(p)
% ------------------------------------------------------------------------------
% 2016-06-12 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A family of CHM obtained from complex Golay sequence;
%
% [1] https://arxiv.org/abs/1204.5164
%
% >> H = G10(rand);
% >> H = G10(0); % ---> BH(10, 4)
% ------------------------------------------------------------------------------

    a = exp(2j * pi * p);
    b = 1/a;

    H = [
             1,      a,    a^2,  i*a^3, -i*a^4,      1,    i*a, -i*a^2,   -a^3,  i*a^4;
        -i*a^4,      1,      a,    a^2,  i*a^3,  i*a^4,      1,    i*a, -i*a^2,   -a^3;
         i*a^3, -i*a^4,      1,      a,    a^2,   -a^3,  i*a^4,      1,    i*a, -i*a^2;
           a^2,  i*a^3, -i*a^4,      1,      a, -i*a^2,   -a^3,  i*a^4,      1,    i*a;
             a,    a^2,  i*a^3, -i*a^4,      1,    i*a, -i*a^2,   -a^3,  i*a^4,      1;
             1, -i*b^4,   -b^3,  i*b^2,   -i*b,     -1, -i*b^4,  i*b^3,   -b^2,     -b;
          -i*b,      1, -i*b^4,   -b^3,  i*b^2,     -b,     -1, -i*b^4,  i*b^3,   -b^2;
         i*b^2,   -i*b,      1, -i*b^4,   -b^3,   -b^2,     -b,     -1, -i*b^4,  i*b^3;
          -b^3,  i*b^2,   -i*b,      1, -i*b^4,  i*b^3,   -b^2,     -b,     -1, -i*b^4;
        -i*b^4,   -b^3,  i*b^2,   -i*b,      1, -i*b^4,  i*b^3,   -b^2,     -b,     -1;
    ];

end
