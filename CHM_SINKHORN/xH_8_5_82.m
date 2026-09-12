function Y = xH_8_5_82(a, b, c, d, f)
% ------------------------------------------------------------------------------
% 2026-09-05 GPT 5.6 Think ==> analytic structure from numerical pattern
%
% five-parameter dephased complex Hadamard matrix of order 8 identified from:
%
% xH_8_5_82_20240201T075439.data
%
% >> H = xH_8_5_L(a, b, c, d, f);
%
% the dependent phase is: e = b + d - a (mod 1)
% ------------------------------------------------------------------------------

    if nargin == 0
        a = 0.15823393394546756;
        b = 0.45875597185198969;
        c = 0.33384398814384963;
        d = 0.031487276475023328;
        f = 0.044070699495149884;
    elseif nargin ~= 5
        error(['xH_8_5_family expects either no arguments or exactly ' ...
               'five arguments: a, b, c, d, f.']);
    endif

    parameters = [a, b, c, d, f];

    if any(!isreal(parameters)) || any(!isfinite(parameters))
        error('All phase parameters must be finite real numbers.');
    endif

%   normalize the five free phases to [0,1).
    a = mod(a, 1);
    b = mod(b, 1);
    c = mod(c, 1);
    d = mod(d, 1);
    f = mod(f, 1);

%   dependent phase
    e = mod(b + d - a, 1);

    A = [
        0  0        0        0        0        0        0        0;
        0  a        b+1/2    a+1/2    c        1/2      b        c+1/2;
        0  d        e        d+1/2    c+1/2    1/2      e+1/2    c;
        0  1/2      1/2      1/2      0        0        1/2      0;
        0  d+1/2    e+1/2    d        c+1/2    1/2      e        c;
        0  f        f+1/2    f        1/2      0        f+1/2    1/2;
        0  f+1/2    f        f+1/2    1/2      0        f        1/2;
        0  a+1/2    b        a        c        1/2      b+1/2    c+1/2;
    ];

    A = mod(A, 1);
    Y = exp(2j*pi*A);

end
