function Y = TU_9_16_3B
% ------------------------------------------------------------------------------
% 2023-01-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% First example from B = BH(9, 3) from Butson Home.
% Diagonal unitary matrix D brings B --> D*B*D' to a 2-unitary form.
% Resulting matrix is in BH(9, 3) again.
%
% Format proof that Y is 2-unitary:
%     Matrix Y is unitary, because after dephasing it reveals the original Butson matrix B.
%     By direct inspection, one can check that dephased forms of Y^R and Y^G
%     are permutation equivalent to B. Hence Y is a 2-unitary matrix.
% ------------------------------------------------------------------------------

    B = B9q3(1);
%%  B = exp(2j*pi*[
%%      0, 0, 0, 0, 0, 0, 0, 0, 0;
%%      0, 0, 0, 1, 1, 1, 2, 2, 2;
%%      0, 0, 0, 2, 2, 2, 1, 1, 1;
%%      0, 1, 2, 0, 1, 2, 0, 1, 2;
%%      0, 1, 2, 1, 2, 0, 2, 0, 1;
%%      0, 1, 2, 2, 0, 1, 1, 2, 0;
%%      0, 2, 1, 0, 2, 1, 0, 2, 1;
%%      0, 2, 1, 1, 0, 2, 2, 1, 0;
%%      0, 2, 1, 2, 1, 0, 1, 0, 2;
%%  ]/3);

    D = diag(exp(2j*pi*[0, 0, 0, 0, 4, 2, 0, 2, 4]/3)); % works also for 1j*pi... and other powers...
    Y = D * B * D';

%    C = exp(2j*pi*[ % = explicit form of Y
%        0, 0, 0, 0, 2, 1, 0, 1, 2;
%        0, 0, 0, 1, 0, 2, 2, 0, 1;
%        0, 0, 0, 2, 1, 0, 1, 2, 0;
%        0, 1, 2, 0, 0, 0, 0, 2, 1;
%        1, 2, 0, 2, 2, 2, 0, 2, 1;
%        2, 0, 1, 1, 1, 1, 0, 2, 1;
%        0, 2, 1, 0, 1, 2, 0, 0, 0;
%        2, 1, 0, 0, 1, 2, 1, 1, 1;
%        1, 0, 2, 0, 1, 2, 2, 2, 2;
%    ]/3);
%    norm(Y - C,"fro");

end

