function H = K9_2(zeta)
% ------------------------------------------------------------------------------
% 2017-03-09 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Matrix K9_2 originally denoted as BC_9^{(2)} - Block Circulant with Circulant Blocks
% taken from the paper of Bengt R. Karlsson;
% [1] https://arxiv.org/pdf/1604.03770.pdf
%
% Watch the scope of (complex) paramater!
% zeta ~= +/-1 is a complex parameter from the intersection of two disks:
% |1 - zeta| <= 4 and |1 + zeta| <= 4
%
% >> H = K9_2(zeta);
% >> H = K9_2(randn + 1j*randn);
% ------------------------------------------------------------------------------

    if ~(abs(1 - zeta) <= 4 && abs(1 + zeta) <= 4) || zeta == 1 || zeta == -1
        warning("No valid ZETA parameter provided! Continue with default zero value!");
        zeta = 0;
    end

    x = 0.25 * (1 + zeta) * (1 + i * sqrt((16 / abs(1 + zeta)^2) - 1));
    y = 0.25 * (1 + zeta) * (1 - i * sqrt((16 / abs(1 + zeta)^2) - 1));
    u = 0.25 * (1 - zeta) * (1 + i * sqrt((16 / abs(1 - zeta)^2) - 1));
    w = 0.25 * (1 - zeta) * (1 - i * sqrt((16 / abs(1 - zeta)^2) - 1));

    H = [
        1 x x  y u w  y w u;
        x 1 x  w y u  u y w;
        x x 1  u w y  w u y;
        
        y w u  1 x x  y u w;
        u y w  x 1 x  w y u;
        w u y  x x 1  u w y;
        
        y u w  y w u  1 x x;
        w y u  u y w  x 1 x;
        u w y  w u y  x x 1;
    ];

end
