function Y = SH_9_16_3
% ------------------------------------------------------------------------------
% 2026-09-10 Claude Opus 5
%
% symmetric complex Hadamard matrix of order 9
%
% source data: 23 independent numerical hits xx_9_16_L_20251213T*.data
% all of which turned out to be monomially EQUIVALENT to each other once dephased correctly
% -- despite d = 16 this is a single isolated matrix
%
% it is a Butson-type matrix BH(9, 3): every core entry is a power of w = exp(2j*pi/3), i.e. a cube root of unity
% this is the catalogued F_3 (x) F_3  matrix (matrix-toolbox.github.io/chm/catalogue/0903B.html)
%
% the log-form L below is the symmetric representative
% ------------------------------------------------------------------------------

    w = exp(2j * pi / 3);
    L = [
        0 0 0 0 0 0 0 0 0;
        0 1 2 0 1 2 0 1 2;
        0 2 1 0 2 1 0 2 1;
        0 0 0 1 1 1 2 2 2;
        0 1 2 1 2 0 2 0 1;
        0 2 1 1 0 2 2 1 0;
        0 0 0 2 2 2 1 1 1;
        0 1 2 2 0 1 1 2 0;
        0 2 1 2 1 0 1 0 2;
    ];

    Y = w .^ L;
end
