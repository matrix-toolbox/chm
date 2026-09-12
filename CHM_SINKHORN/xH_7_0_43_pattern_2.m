function Y = xH_7_0_43_pattern_2
% 2022-09-02
% 2023-01-07

    kMax = 10000;
    do
        printf("RESET\n");
        Z = +Inf;
        Z_OPTIMAL = +Inf;
        mu = 1;

        p1 = rand;
        p3 = rand;

        k = 0;

        while Z_OPTIMAL > 5e-13 && k<kMax
            k++;
            RESTORE_p1 = p1;
            RESTORE_p3 = p3;

            p1 = mod(p1 + randn*mu, 1);
            p3 = mod(p3 + randn*mu, 1);

            a = exp(2j*pi*p1);
            c = exp(2j*pi*p3);

            %        d = a * c / b;
            %        H = [
            %                1   1          1     1        1               1         1     ;
            %                1   a          b     c        b*b/a           b*b       a*b/d ;
            %                1   d          a/d   b        c*b             b*b/d     b*b   ;
            %                1   d*d/a      d/a   c*b*d/a  c               c*b       b*b/a ;
            %                1   d*c        b/a   c*d/a    c*b*d/a         b         c     ;
            %                1   a/b        1/b   b/a      d/a             a/d       b     ;
            %                1   d*c*a/b    a/b   d*c      d*d/a           d         a     ;
            %        ];

            b = (-1 + sqrt(1 - 4*(1 + 1/a + 1/c)*(1 + a + c)) ) / 2 / (1 + 1/a + 1/c);

            Y = [
                1   1             1     1        1            1            1     ;
                1   a             b     c        b*b/a        b*b          b*b/c ;
                1   a*c/b         b/c   b        c*b          b*b*b/a/c    b*b   ;
                1   a*c*c/b/b     c/b   c*c      c            c*b          b*b/a ;
                1   a*c*c/b       b/a   c*c/b    c*c          b            c     ;
                1   a/b           1/b   b/a      c/b          b/c          b     ;
                1   a*a*c*c/b/b   a/b   a*c*c/b  a*c*c/b/b    a*c/b        a     ;
            ];

            Z = nh(Y);
            if Z < Z_OPTIMAL
	        Z_OPTIMAL = Z;
	        mu = Z * 0.01;
                printf("%2.14g \t k=%d\n", Z, k);
                k = 0;
            else
	        p1 = RESTORE_p1;
                p3 = RESTORE_p3;
            end
        end % while Z...
    until (k!=kMax)

    a,c

end
