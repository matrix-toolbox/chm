function Y = Y_9_0_89
% ------------------------------------------------------------------------------
% 2022-12-06 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% An Isolated and symmetric CHM matrix of order 9 with #L(Y) = 89
% found by the Sinkhorn algorithm.
%
% [1] https://arxiv.org/abs/2204.11727
% ------------------------------------------------------------------------------
% 2022-12-06T000337
% Original numerical parameters (included in the general solution)
%    a =     0.7078651580076892 +     0.7063475901273758i;
%    b =     0.0229073020853541 -     0.9997375933269539i;
%    c =    -0.6777758037884053 +     0.7352686310451292i;
%    d =    -0.9393801465461796 +     0.3428774420618481i;
%
% No parameter is a "simple" function of the others.
% ------------------------------------------------------------------------------

    printf("Wait, calculating solution...");
    do
        [u info]=fsolve(@uc, exp(2j*pi*rand(1, 4)));
        for j=1:64
            [u info]=fsolve(@uc, u);
        end
        s = abs(abs(u(1)) - 1) + abs(abs(u(2)) - 1) + abs(abs(u(3)) - 1) + abs(abs(u(4)) - 1);
        printf(".");
    until s < 1e-13
    printf(" Solved!\n");

    a = u(1)
    b = u(2)
    c = u(3)
    d = u(4)

    Y = [
        1 ,   1      ,   1            ,  1           ,  1         ,  1         ,  1         ,   1          ,   1          ;
        1 ,   c^2/b  ,   a            ,  b           ,  c^2       ,  c^2/d     ,  c         ,   c^2/a      ,   d          ;
        1 ,   a      ,   a^2*d^2/c^4  ,  a*b*d/c^3   ,  a*d/c/b   ,  a*d/c^2/b ,  a*b*d/c^2 ,   a*d^2/c^4  ,   d*a/c^2    ;
        1 ,   b      ,   a*b*d/c^3    ,  c*b^2/d     ,  c^3*b/a/d ,  c*b/d     ,  b^2       ,   b*d/c      ,   d/c        ;
        1 ,   c^2    ,   a*d/c/b      ,  c^3*b/a/d   ,  c^3/d     ,  c         ,  c^2*b/d   ,   d/c        ,   d/b        ;
        1 ,   c^2/d  ,   a*d/c^2/b    ,  c*b/d       ,  c         ,  1/c       ,  c^2*b/a/d ,   d/c/b      ,   d/c^2      ;
        1 ,   c      ,   a*b*d/c^2    ,  b^2         ,  c^2*b/d   ,  c^2*b/a/d ,  b^2/c     ,   b          ,   b*d/c^2    ;
        1 ,   c^2/a  ,   a*d^2/c^4    ,  b*d/c       ,  d/c       ,  d/c/b     ,  b         ,   d^2/c^2/b  ,   d^2/c^2    ;
        1 ,   d      ,   d*a/c^2      ,  d/c         ,  d/b       ,  d/c^2     ,  b*d/c^2   ,   d^2/c^2    ,   d/a        ;
    ];

end


function u=uc(v)
    % Mathematica 12.0 has hard times to solve it; even numerically it cannot find all solutions!
    u(1) = 1 + v(1) + v(2) + v(3) + v(3)^2 + v(3)^2/v(1) + v(3)^2/v(2) + v(3)^2/v(4) + v(4);
    u(2) = (v(3) + v(4))/v(3) + (v(2)^2*(v(3) + v(4)))/v(4) +  v(2)*(1 + v(3)/v(4) + v(3)^3/(v(1)*v(4)) + (v(1)*v(4))/v(3)^3 + v(4)/v(3));
    u(3) = 1 + 1/v(1) + (v(2)*v(4))/v(3)^3 + ((1 + 1/v(2) + v(2))*v(4))/v(3)^2 + v(4)/(v(2)*v(3)) + ((1 + v(1))*v(4)^2)/v(3)^4;
    u(4) = 1 + v(4)^2/v(3)^4 + v(4)/v(3) + v(2)/v(3) + v(4)/v(3)^3 + (v(4)^2*v(1))/(v(2)*v(3)^4) + v(2)/v(1) + v(4)^2/(v(2)*v(3)^3) + v(4)/v(3)^2;
end

