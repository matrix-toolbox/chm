function Y = Y_8_0_70
% ------------------------------------------------------------------------------
% 2022-08-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% An isolated matrix found by the Sinkhorn algorithm.
% After permuting it is EQUAL to V_8^0;
% https://chaos.if.uj.edu.pl/~karol/hadamard/catalogue/0804.html
% ------------------------------------------------------------------------------

    a = -0.7799568409992219 + 0.6258333054244768 * I;
    b = -0.6183387264654941 + 0.7859117121871988 * I;
    c =  0.1941704033723143 + 0.9809678151979426 * I;
    d = -0.9991126463144211 + 0.0421179293720193 * I;

    Y = [
        1 ,     1   ,   1    ,   1         ,   1      ,   1    ,   1    ,     1     ;
        1 ,    -1/a ,  -a    ,  -a/b       ,  -b/a    ,   c    ,   1/c  ,    -1     ;
        1 ,    -a   ,   a*a/b,   a*a       ,  -a*c    ,   b    ,   a    ,    -a/c   ;
        1 ,    -a/b ,   a*a  ,   c*a*a/b   ,  -d*c*a/b,   a*c  ,   b/c  ,     a/c   ;
        1 ,    -b/a ,  -a*c  ,  -a*d*c/b   ,  -c*b    ,   d*c  ,  -b/c  ,    -d/a   ;
        1 ,     c   ,   b    ,   a*c       ,   d*c    ,  -b    ,   d    ,     d/a   ;
        1 ,     1/c ,   a    ,   b/c       ,  -b/c    ,   d    ,   d/a/c,     d/c   ;
        1 ,    -1   ,  -a/c  ,   a/c       ,  -d/a    ,   d/a  ,   d/c  ,    -d/c   ;
    ];

end
