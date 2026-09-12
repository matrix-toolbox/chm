function F=FN(N)
% ------------------------------------------------------------------------------
% 2022-12-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Dephased Fourier matrix of size N.
% This script is required by many other scripts: "F*.m", defining particular families stemming from F(N).
%
% >> H = FN(11);
% ------------------------------------------------------------------------------

    F = zeros(N,N);
    for j=0: N-1
    for k=0: N-1
        F(j+1, k+1) = exp(2j*pi*(j+0)*(k+0)/N);
    end, end

end

