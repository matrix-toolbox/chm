function AR = reshuffle(A);
% 20020404 KZ
%
% reshuffling (realignment) process of N = M x M matrices
% technically: for each block from left to right and top to bottom { get block and make it row }
%
% shortcut: R.m

    s1 = size(A); N = s1(1); n2 = s1(2);
    if N ~= n2; NN = N, n22 = n2, else M = fix(sqrt(N));
    if M * M ~= N; mm = M, else AR = A;
    for m = 1 : M;
    for n = 1 : M;
        j = (m - 1) * M + n;
        c = reshape(A(j, :), M, M).';
        ma = (m - 1) * M + 1; mb = m * M;
        na = (n - 1) * M + 1; nb = n * M;
        AR(ma:mb, na:nb) = c;
        end;
        end;
    end;
    end;
end
