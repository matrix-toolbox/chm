function Y = Y_9_0_105
% ------------------------------------------------------------------------------
% 2018-05-06 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Output from the Sinkhorn algorithm:
%    a = -0.3396214473934099 + 0.9405622108454099*1j;
%    b = -0.9635154245261968 + 0.2676528099985182*1j;
%    c = -0.0364845754738033 + 0.9993342162422422*1j;
%    d = +0.8396214473934099 + 0.5431720032153677*1j;
%
% H = [
%    1  1  1  1  1  1  1  1  1  ;
%    1  a  d  a' c' b' c  b  d' ;
%    1  b  c  b' a  d' a' d  c' ;
%    1  c  b' c' d  a  d' a' b  ;
%    1  b  c' b  a' d  a  d' c  ;
%    1  d  a' d' b  c' b' c  a  ;
%    1  a' d' a  c  b  c' b' d  ;
%    1  c' b  c  d' a' d  a  b' ;
%    1  d' a  d  b' c  b  c' a' ;
% ];
%
% 2021-06-09
% Analytic formulas for {a, b, c, d} can be calculated as roots of a monic palindromic polynomial.
% This matrix can be brought to a symmetric form.
%
% [1] https://arxiv.org/abs/2204.11727
% ------------------------------------------------------------------------------

    w = exp(2j * pi * 5 / 6);
    %a=1;
    %z1 = (1/2^((1/3)+a) * w^a * (43+(-1)^a*3*sqrt(771))^(1/3);
    %a=2;
    %z2 = (1/2^((1/3)+a) * w^a * (43+(-1)^a*3*sqrt(771))^(1/3);
    %g1 = sqrt(1-z2-z2')/4 - 1/4;
    %p1 = sqrt(g1^2 - 1) - g1;
    %g2 = sqrt(1-z1-z2')/2 + 1/2;
    %p1 = sqrt(g2^2 - 1) - 2g;

    zeta1 = w/2^(4/3)*(43-3j*sqrt(771))^(1/3);
    zeta2 = w^2*2^(5/3)*(43+3j*sqrt(771))^(1/3);

    gamma1 = sqrt(1-zeta2-conj(zeta2))/4 - 1/4;
    gamma4 = sqrt(1-zeta2-conj(zeta2))/4 + 1/4;
    gamma2 = sqrt(2)*sqrt(1-zeta1-conj(zeta1))/2 + 1/2;
    gamma3 = sqrt(2)*sqrt(1-zeta1-conj(zeta1))/2 - 1/2;

    a = sqrt(gamma1^2 - 1) - gamma1;
    b = sqrt(gamma2^2 - 1) - gamma2;
    c = sqrt(gamma3^2 - 1) + gamma3;
    d = sqrt(gamma4^2 - 1) + gamma4;

    Y = [
        1 1 1 1 1 1 1 1 1;
        1 a d 1/a 1/c 1/b c b 1/d;
        1 b c 1/b a 1/d 1/a d 1/c;
        1 c 1/b 1/c d a 1/d 1/a b;
        1 1/b 1/c b 1/a d a 1/d c;
        1 d 1/a 1/d b 1/c 1/b c a;
        1 1/a 1/d a c b 1/c 1/b d;
        1 1/c b c 1/d 1/a d a 1/b;
        1 1/d a d 1/b c b 1/c 1/a;
    ];

end

