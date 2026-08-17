function Y = Y_10_2_76(p)
% -----------------------------------------------------------------------------
% 2022-08-08 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% A 2-parameter affine family stemming from BH(10, 12) with generic defect = 2 and #L = 76.
% Found quite _accidentally_ using the Sinkhorn algorithm with a random seed.
% It started the series of Y_10_[1,2,3,4]...
%
% [1] https://arxiv.org/abs/2204.11727
% -----------------------------------------------------------------------------

    a = exp(2j * pi * p(1));
    b = exp(2j * pi * p(2));

    w = exp(2j*pi/12);
    Y = [
        1   1   1     1      1      1      1       1       1       1;
        1   1  -i     i      i     -i      w^4     w^4     w^8     w^8;
        1   1   i    -1     -1      i      w^11    w^11    w^7     w^7;
        1   i   w^8   w^10   w^7    w^11   a      -a       w^4     w^4;
        1   i   w^4   w^2    w^11   w^7    w^8     w^8     b      -b;
        1  -i   w^8   w^4    w^7    w^5    a*i    -a*i     w       w;
        1  -i   w^4   w^8    w^11   w      w^5     w^5     b*i    -b*i;
        1  -1   1    -1      1     -1      a*w^11 -a*w^11  b*w^7  -b*w^7;
        1  -1  -i    -i      i      i      a*w^7  -a*w^7   b*w^11 -b*w^11;
        1  -1   i     1     -1     -i     -a       a      -b       b;
    ];


end
