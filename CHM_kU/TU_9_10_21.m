function Y = TU_9_10_21(p)
% ------------------------------------------------------------------------------
% 2022-12-21 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2022-12-24
%
% Two-Unitary 2-parametric family of CHM.
% ------------------------------------------------------------------------------

    w = exp(2j*pi/3);

    a = exp(2j*pi*p(1));
    b = exp(2j*pi*p(2));

    Y = [
        1 ,      1    ,   1    ,   1     ,  1     ,  1     ,  1       , 1       , 1        ;
        1 ,      w    ,   1/w  ,   a/w   ,  a     ,  a/w/w ,  b       , b*w     , b/w      ;
        1 ,      1/w  ,   w    ,   b     ,  b/w   ,  b*w   ,  b*w/a   , b/a     , b*w*w/a  ;
        1 ,      w    ,   1/w  ,   a     ,  a/w/w ,  a/w   ,  b/w     , b       , b*w      ;
        1 ,      1/w  ,   w    ,   b*w   ,  b     ,  b/w   ,  b/a     , b*w*w/a , b*w/a    ;
        1 ,      1    ,   1    ,   w     ,  w     ,  w     ,  1/w     , 1/w     , 1/w      ;
        1 ,      w    ,   1/w  ,   a/w/w ,  a/w   ,  a     ,  b*w     , b/w     , b        ;
        1 ,      1/w  ,   w    ,   b/w   ,  b*w   ,  b     ,  b*w*w/a , b*w/a   , b/a      ;
        1 ,      1    ,   1    ,   1/w   ,  1/w   ,  1/w   ,  w       , w       , w        ;
    ];

    DL = diag([1 1 1 1 a*w b*w w/b w*a/b 1]);
    DR = diag([1 w w w 1 w w w 1]);
    Y = DL * Y * DR;

end

