function Y = Y_10_0_99A
% -----------------------------------------------------------------------------
% 2023-01-07 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% An isolated and symmetric CHM of order 10, with #L = 99
% found by the Sinkhorn algorithm @20220826T092852.
% Parameters (a,b,c,d) are obtained as a solution of the system of non-linear
% equations built out of unitarity constraints.
% Solution:
%
%   a =  0.4894 - 0.8721i, b = -0.5713 - 0.8208i, c = -0.9454 + 0.3259i, d = -0.1074 - 0.9942i
%   a = -0.9454 - 0.3259i, b = -0.9953 + 0.0965i, c = -0.7468 + 0.6650i, d =  0.4255 + 0.9049i
%   a = -0.4085 + 0.9127i, b =  0.9631 - 0.2690i, c =  0.7632 - 0.6462i, d = -0.8146 + 0.5800i
%   :                      :                      :                      :
%
% is in accordance with "random walk" procedure, cf. "Y_10_0_99.m".
%
% [1] https://arxiv.org/abs/2204.11727
% -----------------------------------------------------------------------------

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

    e = a*b*c*d;
    f = -a*(2 + b) - (2*(1 + d + c*d))/d;

    Y = [
	1 , 1   , 1         , 1       , 1       , 1     , 1           , 1     , 1       , 1     ;
	1 , f   , 1/d       , 1/d     , a       , a     , b*a         , c     , c       , 1     ;
	1 , 1/d , 1/d^2/a   , 1/d^2   , 1/d     , a     , e/c^2/d^2   , c^2/e , c/d     , 1/c/d ;
	1 , 1/d , 1/d^2     , 1/c/d^2 , a/c/d   , c*a/e , b/d         , c     , c^2*b/e , 1/d   ;
	1 , a   , 1/d       , a/c/d   , e*a/b/c , a^2   , a^2*b       , a     , e/b     , 1/b   ;
	1 , a   , a         , c*a/e   , a^2     , a^2/c , e*a/c       , a*c   , c       , a/c   ;
	1 , b*a , e/c^2/d^2 , b/d     , a^2*b   , e*a/c , e^2/c^2/d^2 , e     , b*a     , b     ;
	1 , c   , c^2/e     , c       , a       , a*c   , e           , c^2/a , c^2     , c/a   ;
	1 , c   , c/d       , c^2*b/e , e/b     , c     , b*a         , c^2   , c^3*d/e , c*d   ;
	1 , 1   , 1/c/d     , 1/d     , 1/b     , a/c   , b           , c/a   , c*d     , d     ;
    ];

end


function u=uc(v)
    u(1) = 1 + v(1) + (v(3)*v(4))/v(2) + v(1)*v(4)*(2 + 1/v(3) + v(3) + v(4)) + (v(1)^2*v(4)*(v(2) + v(3)*v(4)))/v(3);
    u(2) = 2 + 1/v(2) + v(2) + v(1)/v(3) + v(3)/v(1) + 1/v(4) + 1/(v(3)*v(4)) + v(4) + v(3)*v(4);
    u(3) = v(1)^2*v(2)*v(3)*v(4)*(1 + v(2) + v(4)) + v(3)*(v(2) + v(4) + v(2)*v(4)) + v(1)*v(2)*(1 + v(3)*v(4))^2;
    u(4) = v(3) + v(2)*v(4)*((1 + v(1))*(1 + v(3))*(v(1) + v(3)) + v(1)^2*v(2)*v(3)*v(4));
end

