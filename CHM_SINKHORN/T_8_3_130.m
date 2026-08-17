function Y = T_8_3_130
% ------------------------------------------------------------------------------
% 2022-08-27 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Check "T_8_3_130A.m" for the original "pattern".
%
% It is difficult to express parameters {d,e,f} as functions of {a,b,c}, but
% once {d,e,f} are fixed, it is easy to find such a triplet {a,b,c}
% that Y(a,b,c) is a valid CHM. In particular, for generic values of {d,e,f},
% one obtains: defect(Y) = 3, #L(Y) = 130.
% ------------------------------------------------------------------------------

    kMax = 10000;
    do
        printf("RESET\n");
        Z = +Inf;
        Z_OPTIMAL = +Inf;
        mu = 1;

        P = rand(1, 3);
        k = 0;

        %-----------------------------------------------------------------------
        d = exp(2j*pi*rand);   % FIX ME or left at random
        e = exp(2j*pi*rand);   % FIX ME or left at random
        f = exp(2j*pi*rand);   % FIX ME or left at random
        d = exp(2j*pi*rand);
        e = exp(2j*pi*rand);
        f = exp(2j*pi*rand);
        %-----------------------------------------------------------------------

        while Z_OPTIMAL > 5e-13 && k<kMax
            k++;
            RESTORE_P = P;

            P = mod(P + randn(1,3)*mu, 1);

            a = exp(2j*pi*P(1));
            b = exp(2j*pi*P(2));
            c = exp(2j*pi*P(3));

            p = (a*(b - d)*f)/(a^2 + c*e*f + a*(1 + b + c + d + e + f));
            q = (a*c*d*e*(1 + a + b + c + d + e + 2*f))/((a*d - b*c*e)*f*(1 + a + b + c + d + e + f));
            r = -((a^2*d*e + b*c*d*e*f + a*b*(c*e*f + d*e*f + c*d*(e + f + e*f)) + a*c*d*e*f)/(a*c*(d - e)*f));
            s = -(a/(1 + a + b + c + d + e + f));
            t = -((1 + d - p - q + r + s)/(1 + r));

            Y=[
                1  ,  1    ,   1    ,   1    ,     1     ,     1    ,     1   ,     1    ;
                1  ,  a    ,   b    ,   c    ,     d     ,     e    ,     f   ,     a/s  ;
                1  ,  a/b  ,  -1    ,   q    ,    -q     ,    -e    ,    -a/b ,     e    ;
                1  ,  a/c  ,   p    ,  -1    ,    -p     ,    -a/c  ,    -f   ,     f    ;
                1  ,  a/d  ,   a/e  ,   a/f  ,     s     ,     a/c  ,     a/b ,     a    ;
                1  ,  a/e  ,  -a/e  ,  -c    ,     t     ,    -1    ,    -t   ,     c    ;
                1  ,  a/f  ,  -b    ,  -a/f  ,     r     ,    -r    ,    -1   ,     b    ;
                1  ,  s    ,  -p    ,  -q    ,     r*t   ,     r    ,     t   ,     d    ;
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

    a, b, c

end

%
% a = -3/5 + 4j/5;
% b = -1;
% c = -a;
%
