function Y=T_8_3_74(p)
% ------------------------------------------------------------------------------
% 2022-08-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 3-parameter family obtained from a random seed applied to the Sinkhorn algorithm.
% Generic values: d(Y) = 3, #L(Y) = 74.
% It stems from BH(8, 2).
%
% Watch out of the scope of parameter;
% Example of call:
% >> Y = T_8_3_74([3/5 1/3 4/9]);
%
% [1] https://arxiv.org/abs/2204.11727
% ------------------------------------------------------------------------------


    a = exp(2j*pi * p(1));
    b = exp(2j*pi * p(2));
    c = exp(2j*pi * p(3));

    if abs(a + b + c + a * b * c) < 1e-10
        warning("Singularity! Continue with (a,b,c)=(1,1,1)...");
        a = 1;
        b = 1;
        c = 1;
    end
    d = (1 + b*c + a*b + a*c) / (a + b + c + a * b * c);

    Y = [
        1  ,      1    ,   1    ,    1   ,    1    ,   1   ,  1       ,    1       ;
        1  ,     -1    ,   c    ,   -c   ,    d    ,  -d   ,  c*d     ,   -c*d     ;
        1  ,      a    ,  -a    ,   -1   ,    a*d  ,   d   , -d       ,   -a*d     ;
        1  ,     -a    ,  -a*c  ,    c   ,    a    ,  -1   , -c       ,    a*c     ;
        1  ,      b    ,  -c    ,   -b*c ,   -1    ,  -b   ,  c       ,    b*c     ;
        1  ,     -b    ,  -1    ,    b   ,   -d    ,   b*d ,  d       ,   -b*d     ;
        1  ,      a*b  ,   a*c  ,    b*c ,   -a*d  ,  -b*d , -c*d     ,   -a*b*c*d ;
        1  ,     -a*b  ,   a    ,   -b   ,   -a    ,   b   , -1       ,    a*b     ;
    ];

end
