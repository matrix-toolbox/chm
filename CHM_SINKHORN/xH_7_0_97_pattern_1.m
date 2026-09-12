function Y = xH_7_0_97_pattern_1
% 2022-09-03
% 2023-01-07

    kMax = 10000;
    do
        printf("RESET\n");
        Z = +Inf;
        Z_OPTIMAL = +Inf;
        mu = 1;

        p = rand(1,6);
        k = 0;

        while Z_OPTIMAL > 5e-13 && k<kMax
            k++;
            RESTORE_p = p;

            p = mod(p + randn(1,6)*mu, 1);

            a = exp(2j*pi*p(1));
            b = exp(2j*pi*p(2));
            c = exp(2j*pi*p(3));
            d = exp(2j*pi*p(4));
            e = exp(2j*pi*p(5));
            f = exp(2j*pi*p(6));

            Y = [
                1    1    1    1    1    1    1;
                1    a    b    c    d    e    f;
                1    c/a  d/a  b/a  e/a  1/a  f/a;
                1    d/b  1/b  e/b  a/b  c/b  f/b;
                1    b/c  e/c  d/c  1/c  a/c  f/c;
                1    e/d  a/d  1/d  c/d  b/d  f/d;
                1    1/e  c/e  a/e  b/e  d/e  f/e;
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

    a,b,c,d,e,f

end
