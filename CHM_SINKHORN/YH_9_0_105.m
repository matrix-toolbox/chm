function Y = YH_9_0_105
% ------------------------------------------------------------------------------
% 2026-09-10 Claude Opus 5
%
% isolated complex Hadamard matrix of order 9
%
% source data: x_9_0_L_20251213T181223.data
%
% identified as Y9^(0) [matrix-toolbox.github.io/chm/catalogue/0906.html]
%
% a, b, c, d are unimodular numbers built from two cubic irrationalities zeta1, zeta2
% a*, b*, c*, d* denote their complex conjugates ( here |a|=|b|=|c|=|d|=1 so a* = 1/a and so on)
% ------------------------------------------------------------------------------

    omega  = exp(1j * pi * 5 / 3);
    s771   = sqrt(771);

    zeta1 = omega / 2^(4/3) * (43 - 3j*s771)^(1/3);
    zeta2 = 2^(5/3) * omega^2  * (43 + 3j*s771)^(1/3);

    gamma1 = (1/4)       * sqrt(1 - zeta2 - conj(zeta2)) - 1/4;
    gamma2 = (sqrt(2)/2) * sqrt(1 - zeta1 - conj(zeta1)) + 1/2;
    gamma3 = (sqrt(2)/2) * sqrt(1 - zeta1 - conj(zeta1)) - 1/2;
    gamma4 = (1/4)       * sqrt(1 - zeta2 - conj(zeta2)) + 1/4;

    a = sqrt(gamma1^2 - 1) - gamma1;
    b = sqrt(gamma2^2 - 1) - gamma2;
    c = sqrt(gamma3^2 - 1) + gamma3;
    d = sqrt(gamma4^2 - 1) + gamma4;

    Y = [
        1  1   1   1   1   1   1   1   1  ;
        1  a   d   a'  c'  b'  c   b   d' ;
        1  b   c   b'  a   d'  a'  d   c' ;
        1  c   b'  c'  d   a   d'  a'  b  ;
        1  b'  c'  b   a'  d   a   d'  c  ;
        1  d   a'  d'  b   c'  b'  c   a  ;
        1  a'  d'  a   c   b   c'  b'  d  ;
        1  c'  b   c   d'  a'  d   a   b' ;
        1  d'  a   d   b'  c   b   c'  a' ;
    ];

end
