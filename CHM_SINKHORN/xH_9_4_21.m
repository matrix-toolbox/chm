function Y = xH_9_4_21(a, b)
% ------------------------------------------------------------------------------
% 2026-09-10 Claude Opus 5
%
% two-parameter dephased complex Hadamard matrix of order 9 identified from xx_9_4_L_20251213T183218.data
%
% >> H = xH_9_4_21(a, b);
%
% the core splits into a:
%     fixed 5x5 Butson(_, 3) block (rows 1-5, entries  multiples of 1/3 only)
%     and a 3x8 free block (rows 6-8)
% where the two free real phases a, b each appear cyclically shifted by 0, +1/3, +2/3 across the three free rows
% i.e. multiplied, row by row, by the three cube roots of  unity in cyclic order
% this is the signature of a Dita-type (generalised tensor) deformation of a Butson(9, 3) skeleton by two extra free phases
% ------------------------------------------------------------------------------

    if nargin == 0
        a = 0.06625401369684453;
        b = 0.164323384108616;
    elseif nargin ~= 2
        error('xH_9_4_21 expects either no arguments or exactly two arguments: a, b!');
    endif

    if any(!isreal([a, b])) || any(!isfinite([a, b]))
        error('all phase parameters must be finite real numbers!');
    endif

    a = mod(a, 1);
    b = mod(b, 1);

    A = [
        2/3,   1/3,   2/3,   2/3,   1/3,   0,     1/3,   0   ;
        1/3,   2/3,   1/3,   1/3,   2/3,   0,     2/3,   0   ;
        0,     0,     1/3,   2/3,   2/3,   2/3,   1/3,   1/3 ;
        2/3,   1/3,   0,     1/3,   0,     2/3,   2/3,   1/3 ;
        1/3,   2/3,   2/3,   0,     1/3,   2/3,   0,     1/3 ;
        a,     b,     a+2/3, a+1/3, b+1/3, 1/3,   b+2/3, 2/3 ;
        a+2/3, b+1/3, a+1/3, a,     b+2/3, 1/3,   b,     2/3 ;
        a+1/3, b+2/3, a,     a+2/3, b,     1/3,   b+1/3, 2/3 ;
    ];

    A = mod(A, 1);

    Y = [ones(1, 9); [ones(8, 1), exp(2j * pi * A)]];

end
