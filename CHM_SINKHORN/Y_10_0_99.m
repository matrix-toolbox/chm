function Y = Y_10_0_99
% -----------------------------------------------------------------------------
% 2023-01-06 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% An isolated and symmetric CHM of order 10, with #L = 99
% found by the Sinkhorn algorithm @20220826T092852.
% See also "Y_10_9_99A.m".
%
% [1] https://arxiv.org/abs/2204.11727
%
% Random walk converges to different solutions, e.g.:
% a = -0.9454 - 0.3259i, b = -0.9953 + 0.0964i, c = -0.7468 + 0.6650i, d =  0.4255 + 0.9049i
% a =  0.4894 + 0.8721i, b = -0.5713 + 0.8208i, c = -0.9454 - 0.3259i, d = -0.1074 + 0.9942i
% a = -0.7468 - 0.6650i, b =  0.9725 + 0.2331i, c =  0.4894 + 0.8721i, d = -0.9196 - 0.3929i
% a = -0.9016 - 0.4326i, b = -0.7520 + 0.6592i, c = -0.4085 - 0.9127i, d = -0.1966 - 0.9805i
% a = -0.4085 - 0.9127i, b =  0.9631 + 0.2690i, c =  0.7632 + 0.6462i, d = -0.8146 - 0.5800i
% a =  0.7632 + 0.6462i, b = -0.1480 + 0.9890i, c = -0.9016 + 0.4326i, d = -0.2469 + 0.9690i
% :                      :                      :                      :
% -----------------------------------------------------------------------------

    kMax = 10000;
    nP = 4;

    do
        printf("RESET\n");
        Z = +Inf;
        Z_OPTIMAL = +Inf;
        mu = 1;
        P = rand(1, nP);

        k = 0;
        while Z_OPTIMAL > 5e-13 && k<kMax
            k++;
            RESTORE_P = P;
            P = mod(P + randn(1, nP)*mu, 1);

            a = exp(2j*pi*P(1));
            b = exp(2j*pi*P(2));
            c = exp(2j*pi*P(3));
            d = exp(2j*pi*P(4));
            e = a*b*c*d;
            f = -a*(2 + b) - 2*(1+d+c*d)/d;

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

            Z = nh(Y);

            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z * 0.01;
                printf("%2.14g \t k=%d\n", Z, k);
                k = 0;
            else
                P = RESTORE_P;
            end
        end % while Z...
    until (k!=kMax)

    a,b,c,d

end

