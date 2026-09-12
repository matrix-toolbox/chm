function Y = xH_10_2_L(p1, p2)
% ------------------------------------------------------------------------------
% 2026-09-10 Claude Opus 5
%
% analytic form derived from numerical data
%
% two-parameter affine family of complex Hadamard matrices of order 10
%
% >> H = xH_10_2(p1, p2); % a = exp(2j * pi * p1) and b = exp(2j * pi * p2)
%
% everything else is a fixed 12th root of unity w = exp(2j * pi / 12)
% parameters a, b are free anywhere on the unit circle
% ------------------------------------------------------------------------------
% p1    p2     d   #L   q
% 0     0      3   12   12
% *     *      2*  76   -
% 0     *      3   36   -
% *     0      3   36   -
% 1/4   1/4    4   12   12
% 1/5   1/6    2   36   60
% 1/5   1/5    2   40   60
% 1/7   1/7    2   40   84
% 1/8   1/8    2   24   24
% 1/9   1/9    2   36   36
% ...
% 1/11  1/121  2   76   1452
% ...
% ------------------------------------------------------------------------------

    if nargin == 0
        p1 = rand();
        p2 = rand();
    endif

    a = exp(2j * pi * p1);
    b = exp(2j * pi * p2);
    w = exp(2j * pi / 12);
    i = 1j;

    Y = [
        1,   1,      1,      1,      1,      1,      1,        1,        1,        1        ;
        1,   1,     -i,      i,      i,     -i,      w^4,      w^4,      w^8,      w^8      ;
        1,   1,      i,     -1,     -1,      i,      w^11,     w^11,     w^7,      w^7      ;
        1,   i,      w^8,    w^10,   w^7,    w^11,   a,       -a,        w^4,      w^4      ;
        1,   i,      w^4,    w^2,    w^11,   w^7,    w^8,      w^8,      b,       -b        ;
        1,  -i,      w^8,    w^4,    w^7,    w^5,    a*i,     -a*i,      w,        w        ;
        1,  -i,      w^4,    w^8,    w^11,   w,      w^5,      w^5,      b*i,     -b*i      ;
        1,  -1,      1,     -1,      1,     -1,      a*w^11,  -a*w^11,   b*w^7,   -b*w^7    ;
        1,  -1,     -i,     -i,      i,      i,      a*w^7,   -a*w^7,    b*w^11,  -b*w^11   ;
        1,  -1,      i,      1,     -1,     -i,     -a,        a,       -b,        b        ;
    ];

end
