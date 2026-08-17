function Y=Y_8_5_42(p)
% ------------------------------------------------------------------------------
% 2022-12-10 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 3-parametric output from the Sinkhorn algorithm
% with generic values: d = 5, #L = 42.
% It stems from BH(8, 2) and it is very likely that it is a part of F8.
% ------------------------------------------------------------------------------

    %a = 0.894356497722562 + 0.447354954126402i;
    %b = 0.132554945206970 + 0.991175658751352i;
    %c = 0.892248459724733 - 0.451544777534680i;

    a = exp(2j*pi * p(1));
    b = exp(2j*pi * p(2));
    c = exp(2j*pi * p(3));

    d = -a*b; % 0.324855957733842 - 0.945763504648403i;
    e = -c/b; % 0.329288244526594 + 0.944229448818769i;

    Y = [
       1 ,  1 ,  1  ,   1 ,  1 ,  1  ,   1  ,    1    ;
       1 ,  1 ,  a  ,   e , -e , -a  ,  -1  ,   -1    ;
       1 ,  1 , -1  ,  -e ,  e , -1  ,   a  ,   -a    ;
       1 ,  1 , -a  ,  -1 , -1 ,  a  ,  -a  ,    a    ;
       1 , -1 ,  b  ,   c , -c ,  b  ,  -b  ,   -b    ;
       1 , -1 , -b  ,   b ,  b , -b  ,   d  ,   -d    ;
       1 , -1 ,  d  ,  -c ,  c , -d  ,  -d  ,    d    ;
       1 , -1 , -d  ,  -b , -b ,  d  ,   b  ,    b    ;
    ];

end
