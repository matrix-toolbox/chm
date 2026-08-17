function Y = T_8_3_130A
% ------------------------------------------------------------------------------
% 2022-08-23 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2023-01-07
%
% Generic values: defect = 3, #L = 130.
% This is an alternative version of "Y_8_3_130.m".
% ------------------------------------------------------------------------------

    kMax = 10000;
    do
        printf("RESET\n");
        Z = +Inf;
        Z_OPTIMAL = +Inf;
        mu = 1;

        p = rand(1,3);
        k = 0;

        %-----------------------------------------------------------------------
        c_1 = exp(2j*pi*rand);
        c_2 = exp(2j*pi*rand);
        c_3 = exp(2j*pi*rand);
        %-----------------------------------------------------------------------

        while Z_OPTIMAL > 5e-13 && k<kMax
            k++;
            RESTORE_p = p;

            p = mod(p + randn(1,3)*mu, 1);

            a = exp(2j*pi*p(1));
            b = exp(2j*pi*p(2));
            c = exp(2j*pi*p(3));

            d = -(c_3 + a + c + b + 1 + c_2 + a/c_1);
            e =  (c_3/c_2 - c_1/c) / (1/d - 1/b);
            f =  (a/c_2 - a/c) / (c_1/c_3 - d/b);
            g = -(1 + a/c_2 + b + f + e) / (e/f + e + 1);

            Y = [
                1 , 1     , 1     ,  1    , 1      ,  1     ,  1     ,  1    ;
                1 , c_1   , c_3   , -c_3  , e      , -1     , -c_1   , -e    ;
                1 , a     , c_2   ,  c_3  , b      ,  a/c_1 ,  c     ,  d    ;
                1 , a/b   , a     ,  a/d  , a/c_2  ,  a/c_3 ,  c_1   ,  a/c  ;
                1 , a/c_2 , b     ,  f    , e*g    ,  g     ,  e*g/f ,  e    ;
                1 , a/d   , c     , -a/d  , g      , -g     , -c     , -1    ;
                1 , a/c   , a/c_1 , -f    , f      , -a/c_1 , -1     , -a/c  ;
                1 , a/c_3 , d     , -1    , e*g/f  , -a/c_3 , -e*g/f , -d    ;
            ];
            Z = nh(Y);
            if Z < Z_OPTIMAL
	        Z_OPTIMAL = Z;
	        mu = Z * 0.01;
                printf("%2.14g \t k=%d\n", Z, k);
                k = 0;
            else
	        p = RESTORE_p;
            end
        end % while Z...
    until (k!=kMax)

    c_1, c_2, c_3, a, b, c

end
