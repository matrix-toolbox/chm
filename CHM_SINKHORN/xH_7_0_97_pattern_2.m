function Y = xH_7_0_97_pattern_2
% 2022-09-03
% 2023-01-07
%
% Symmetric pattern. Recovers symmetric form of Q7.

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
                1	1	1	1	1	1	1	;
                1	a	b	c	d	e	f	;
                1	b	c	e	f	d	a	;
                1	c	e	d	a	f	b	;
                1	d	f	a	c	b	e	;
                1	e	d	f	b	a	c	;
                1	f	a	b	e	c	d	;
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


%
% Clear[a];
% Clear[b];
% Clear[c];
% Clear[d];
% Clear[e];
% Clear[f];
% NSolve[{
%   1 +   a +   b +   c +   d +   e +   f == 0,
%   1 + b/a + c/b + e/c + f/d + d/e + a/f == 0,
%   1 + c/a + e/b + d/c + a/d + f/e + b/f == 0,
%   1 + d/a + f/b + a/c + c/d + b/e + e/f == 0,
%   1 + e/a + d/b + f/c + b/d + a/e + c/f == 0,
%   1 + f/a + a/b + b/c + e/d + c/e + d/f == 0
%  }, {a, b, c, d, e, f}]
