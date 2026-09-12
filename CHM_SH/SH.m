function Y = SH(X, km)
% ------------------------------------------------------------------------------
% 2023-01-14
% modification of Sinkhorn algorithm for symmetric CHM
% ------------------------------------------------------------------------------

    error('OBSOLETE!')

    N = size(X, 1);
    k = 0;
    while (nh(X) >= 1e-13 || n1(X) >= 5e-15) && (k < km)
        X = X ./ abs(X);
        X = PD_SVDi(X);
        X = (X + X') / 2;
        X = X * sqrt(N);
        k++;
     end;

    Y = X; % |H| = 1, H'*H = N*I
end

#            fileName = strcat('X', int2str(N), '_', int2str(ud(X,'S',1e-8)), '_', int2str(L), '_', datestr(now(), 30), '.dat');
#            CHMS = struct('X0', X0, 'H', H, 'K', K, 'L', L, 'defect', udefect, 'k', k);
#            save(fileName, 'CHMS')

