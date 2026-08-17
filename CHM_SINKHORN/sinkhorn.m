function Y=sinkhorn(X, kMax)
% 2022-04-06
% 2023-01-14
% Modification of the Sinkhorn's algorithm:
% consecutive projections onto unimodular (|H_jk=1|) and orthogonal (H*H'=Id_N) subspaces of C(NxN).
% Once successfully used by V. Elser (2011) for N=8, which resulted with a new class of CHM.

    addpath ../matrix_tool

    N = size(X, 1);
    k = 0;
    while (nh(X)>=1e-13 || n1(X)>=5e-15) && (k<kMax)
        X = X ./ abs(X);
        X = PD(X);
        X = X * sqrt(N);
        k++;
     end;

    Y = X; % |H| = 1, H'*H = N*I
end




#            fileName = strcat('X', int2str(N), '_', int2str(ud(X,'S',1e-8)), '_', int2str(L), '_', datestr(now(), 30), '.dat');
#            CHMS = struct('X0', X0, 'H', H, 'K', K, 'L', L, 'defect', udefect, 'k', k);
#            save(fileName, 'CHMS')

