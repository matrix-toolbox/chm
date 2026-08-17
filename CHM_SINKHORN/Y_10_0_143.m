function Y = Y_10_0_143
% -----------------------------------------------------------------------------
% 2023-01-09 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% An isolated CHM found by the Sinkhorn algorithm.
%
% [1] https://arxiv.org/abs/2204.11727
%
% General solution contains D_10, N_10 and probably some other matrices.
% At least one particular solution provides an isolated CHM
% of order N = 10 with #L = 143 (not BH!)
% The convergence to this particular one is not so frequent...
% ------------------------------------------------------------------------------
% a = -1
% b =  1
% c = -1
% d =  1
% e = -I
% f = -1
% g =  I
% h =  I ==> d = 16, #L = 4, BH(10, 4)
% ------------------------------------------------------------------------------
% a = -1
% b =  I
% c = -I
% d =  1
% e = -1
% f = -I
% g =  I
% h =  I ==> d = 11, #L = 4, BH(10, 4)
% ------------------------------------------------------------------------------
% a = -0.8750 - 0.4841i
% b = -0.8750 - 0.4841i
% c =  1
% d =  1
% e = -0.8750 + 0.4841i
% f = -0.2500 - 0.9682i
% g = -0.2500 + 0.9682i
% h = -0.2500 - 0.9682i ==> d = 0, #L = 9, not BH (Nicoara)
% ------------------------------------------------------------------------------
% a = -0.8889 - 0.4580i
% b = -0.7206 + 0.6933i
% c =  0.2128 - 0.9771i
% d =  0.7052 + 0.7090i
% e = -0.9997 - 0.0220i
% f =  0.3590 - 0.9333i
% g = -0.2491 + 0.9685i
% h =  0.6182 + 0.7860i ==> d = 0, #L = 143 <------------------- this is the one
% ------------------------------------------------------------------------------

    kMax = 10000;
    nP = 8;

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
            e = exp(2j*pi*P(5));
            f = exp(2j*pi*P(6));
            g = exp(2j*pi*P(7)); % simple function of {a}
            h = exp(2j*pi*P(8)); % simple function of {a, d, g}

            Y = [
                1 , 1    , 1    , 1  , 1  , 1      , 1  , 1  , 1    , 1     ;
                1 , 1    , 1    , 1  , a/g, g/a    , a  , a  , 1/a  , 1/a   ;
                1 , b/a^2, c/a^2, d/a, 1/a, g/a    , b/a, c/a, e/a  , f/a   ;
                1 , b/a^2, c/a^2, d/a, 1  , g/a^2  , e  , f  , b/a^2, c/a^2 ;
                1 , g/a^2, g/a^2, g/a, g/a, g^2/a^2, g  , g  , g/a  , g/a   ;
                1 , e    , f    , d  , a  , g      , b  , c  , b/a  , c/a   ;
                1 , f    , e    , d  , a  , g      , c  , b  , c/a  , b/a   ;
                1 , d/a  , d/a  , h  , 1  , g/a    , d  , d  , d/a  , d/a   ;
                1 , c/a^2, b/a^2, d/a, 1  , g/a^2  , f  , e  , c/a^2, b/a^2 ;
                1 , c/a^2, b/a^2, d/a, 1/a, g/a    , c/a, b/a, f/a  , e/a   ;
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

    a,b,c,d,e,f,g,h

end

