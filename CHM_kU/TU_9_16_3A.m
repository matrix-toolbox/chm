function Y = TU_9_16_3A
% ------------------------------------------------------------------------------
% 2022-12-29 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Original diagonals for F3 (x) F3:
%	DL = diag(exp(2j*pi*[0, 1/3, 1/3, 1/3, 1/3, 0, 0, 2/3, 0]));
%	DR = diag(exp(2j*pi*[0, 1/3, 1/3, 1/3, 0, 1/3, 2/3, 2/3, 1/3]));
% ------------------------------------------------------------------------------

    Y = kron(F3, F3);

%   diagonals that work for F3(x)F3 and also for TU_9_4_57:
    DL = diag(exp(2j*pi*[0, 1, 1, 1, 1, 0, 0, 2, 0]/3));  
    DR = diag(exp(2j*pi*[0, 1, 1, 1, 0, 1, 2, 2, 1]/3));

    Y = DL * Y *DR;

end

