function H = B9
% ------------------------------------------------------------------------------
% 2023-02-04 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated CHM of order N = 9 with defect = 2, found by K. Beauchamp and R. Nicoara;
% [1] https://arxiv.org/abs/math/0609076
% [2] https://arxiv.org/abs/1604.03770
% [3] https://arxiv.org/abs/1609.02829
%
% >> H = B9;
% ------------------------------------------------------------------------------

    H = exp(2j*pi*[
        0 0 0 0 0 0 0 0 0;
        0 7 9 4 9 5 5 1 3;
        0 9 1 5 5 3 2 7 7;
        0 4 5 7 1 3 9 9 5;
        0 9 5 1 3 5 7 6 1;
        0 5 3 3 5 9 7 1 8;
        0 5 2 9 7 7 1 5 3;
        0 1 7 9 6 1 5 3 5;
        0 3 7 5 1 8 3 5 9;
    ]/10);

%    ORIGINAL non-symmetric version:
%
%    e = exp(2 * pi * i / 10);
%
%    H = [
%        1,   1,   1,   1,   1,   1,   1,   1,   1;
%        1,  -1, e^3, e^3,  -1, e^9, e^8, e^7,   e;
%        1, e^4,  -1, e^7,   e, e^3,  -1, e^9, e^9;
%        1, e^3, e^7,  -1,   e, e^8, e^9, e^3,  -1;
%        1, e^9,   e,  -1,  -1, e^3, e^7, e^2, e^7;
%        1, e^9,  -1,   e, e^3,  -1,   e, e^7, e^6;
%        1,   e, e^7, e^9, e^6,   e,  -1,  -1, e^3;
%        1, e^7, e^9, e^4, e^9,  -1, e^3,  -1,   e;
%        1,  -1, e^2, e^9, e^7, e^7, e^3,   e,  -1;
%    ];

end
