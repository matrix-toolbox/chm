function H=MH36
% ------------------------------------------------------------------------------
% 2022-07-10 Padraig O' Cathain
% 2023-02-04 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Menon-Hadamard with cyclic property.
%
% SL = linear entropy
% R = reshuffling
% T = partial transpose
%
% >> SL(H) = 2/3 = 0.(6)
% >> SL(H^R) = 1
% >> SL(H^T) = 1 ---> dual unitary matrix
% ------------------------------------------------------------------------------

    n = -1;
    H=[
    1 n n n n n,n 1 n n n n,n n 1 n n n,n n n 1 n n,n n n n 1 n,n n n n n 1;
    n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n,n 1 1 1 1 n;
    n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n,n 1 1 1 n 1,1 n 1 1 1 n;
    n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n,n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n;
    n 1 1 1 n 1,1 n 1 1 1 n,n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n;
    n 1 1 1 1 n,n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n;

    n 1 1 1 1 n,n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n;
    1 n n n n n,n 1 n n n n,n n 1 n n n,n n n 1 n n,n n n n 1 n,n n n n n 1;
    n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n,n 1 1 1 1 n;
    n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n,n 1 1 1 n 1,1 n 1 1 1 n;
    n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n,n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n;
    n 1 1 1 n 1,1 n 1 1 1 n,n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n;

    n 1 1 1 n 1,1 n 1 1 1 n,n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n;
    n 1 1 1 1 n,n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n;
    1 n n n n n,n 1 n n n n,n n 1 n n n,n n n 1 n n,n n n n 1 n,n n n n n 1;
    n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n,n 1 1 1 1 n;
    n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n,n 1 1 1 n 1,1 n 1 1 1 n;
    n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n,n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n;

    n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n,n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n;
    n 1 1 1 n 1,1 n 1 1 1 n,n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n;
    n 1 1 1 1 n,n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n;
    1 n n n n n,n 1 n n n n,n n 1 n n n,n n n 1 n n,n n n n 1 n,n n n n n 1;
    n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n,n 1 1 1 1 n;
    n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n,n 1 1 1 n 1,1 n 1 1 1 n;

    n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n,n 1 1 1 n 1,1 n 1 1 1 n;
    n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n,n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n;
    n 1 1 1 n 1,1 n 1 1 1 n,n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n;
    n 1 1 1 1 n,n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n;
    1 n n n n n,n 1 n n n n,n n 1 n n n,n n n 1 n n,n n n n 1 n,n n n n n 1;
    n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n,n 1 1 1 1 n;

    n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n,n 1 1 1 1 n;
    n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n,n 1 1 1 n 1,1 n 1 1 1 n;
    n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n,n 1 1 n 1 1,1 n 1 1 n 1,1 1 n 1 1 n;
    n 1 1 1 n 1,1 n 1 1 1 n,n 1 n 1 1 1,1 n 1 n 1 1,1 1 n 1 n 1,1 1 1 n 1 n;
    n 1 1 1 1 n,n n 1 1 1 1,1 n n 1 1 1,1 1 n n 1 1,1 1 1 n n 1,1 1 1 1 n n;
    1 n n n n n,n 1 n n n n,n n 1 n n n,n n n 1 n n,n n n n 1 n,n n n n n 1;
    ];

end
