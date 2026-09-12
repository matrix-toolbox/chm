function Y = xH_7_0_43_pattern_1
% 2022-09-03
% 2023-01-07

    kMax = 10000;
    do
        printf("RESET\n");
        Z = +Inf;
        Z_OPTIMAL = +Inf;
        mu = 1;

        p1 = rand;
        p2 = rand;
        p3 = rand;

        k = 0;

        while Z_OPTIMAL > 5e-13 && k<kMax
            k++;
            RESTORE_p1 = p1;
            RESTORE_p2 = p2;
            RESTORE_p3 = p3;

            p1 = mod(p1 + randn*mu, 1);
            p2 = mod(p2 + randn*mu, 1);
            p3 = mod(p3 + randn*mu, 1);

            a = exp(2j*pi*p1);
            b = exp(2j*pi*p2);
            c = exp(2j*pi*p3);

            %        d = b * c / a;
            %        Y = [
            %                1       1       1        1      1      1       1      ;
            %                1       a       b        c      a*a/c  a*a     a*a/b  ;
            %                1       c       a*a/c    a      b      a*a/b   a*a    ;
            %                1       d       b*b      d*d/c  b*b/d  b       a*a/c  ;
            %                1       d*d     d        d*d/a  d*d/c  a       c      ;
            %                1       d*d/c   b*b/d    d      b*b    a*a/c   b      ;
            %                1       d*d/a   d*d/c    d*d    d      c       a      ;
            %        ];

            Y = [
                1       1               1           1               1           1       1       ;
                1       a               b           c               a*a/c       a*a     a*a/b   ;
                1       c               a*a/c       a               b           a*a/b   a*a     ;
                1       b*c/a           b*b         b*b*c/a/a       a*b/c       b       a*a/c   ;
                1       b*b*c*c/a/a     b*c/a       b*b*c*c/a/a/a   b*b*c/a/a   a       c       ;
                1       b*b*c/a/a       a*b/c       b*c/a           b*b         a*a/c   b       ;
                1       b*b*c*c/a/a/a   b*b*c/a/a   b*b*c*c/a/a     b*c/a       c       a       ;
            ];

            Z = nh(Y);
            if Z < Z_OPTIMAL
	        Z_OPTIMAL = Z;
	        mu = Z * 0.01;
                printf("%2.14g \t k=%d\n", Z, k);
                k = 0;
            else
	        p1 = RESTORE_p1;
	        p2 = RESTORE_p2;
                p3 = RESTORE_p3;
            end
        end % while Z...
    until (k!=kMax)

    a,b,c

end
