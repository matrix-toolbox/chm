function H = N10B(p)
% ------------------------------------------------------------------------------
% 2016-06-12 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Extension of 1-parameter family found by K. Beauchamp and R. Nicoara
% to a 3-parameter family by P. Lampio et al.;
% [1] https://arxiv.org/abs/1204.5164
%
% >> H = N10B(rand(1, 3));
% ------------------------------------------------------------------------------
    a = exp(2j * pi * p(1));
    b = exp(2j * pi * p(2));
    c = exp(2j * pi * p(3));

    H = [
        1,  1  ,  1      ,  1      ,  1  ,  1,  1,        1,        1,    1;
        1,  1  , -1      , -1      ,  1  ,  1, -i,       -1,       -1,    i;
        1,  a  ,  c      , -i*c    , -a  , -1, -c, -i*a*b*c,  i*a*b*c,  i*c;
        1, -a  , -i*c    ,  c      ,  a  , -1, -c,  i*a*b*c, -i*a*b*c,  i*c;
        1, -i  , -i*c/a  ,  i*c/a  , -1  ,  1,  i,      b*c,     -b*c,   -1;
        1,  i/b, -i*c/a/b,  i*c/a/b, -i/b, -1,  c,       -c,      i*c, -i*c;
        1, -i/b,  i*c/a/b, -i*c/a/b,  i/b, -1,  c,      i*c,       -c, -i*c;
        1, -1  ,  i*c/a  , -i*c/a  , -i  ,  1,  i,     -b*c,      b*c,   -1;
        1,  i  , -c      , -c      ,  i  , -i, -1,        c,        c,   -i;
        1, -1  ,  i*c    ,  i*c    , -1  ,  i, -i,     -i*c,     -i*c,    1;
    ];

end
