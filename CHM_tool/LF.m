function L = LF(X, BHF_detector = 0)
% ------------------------------------------------------------------------------
% 2023-02-03
%
% If X in BH(N, q) then print LOG-form of X.
% ------------------------------------------------------------------------------
% 2024-06-16
%
% 1.9999999999994 ---> 2
% 2.0000000000004 ---> 2
% q-1.99999999999 ---> 0
% ------------------------------------------------------------------------------

    N = size(X, 1);
    [q _] = isBH(X, 10000, '');
    L = zeros(N, N);
    if q
        A = mod(round(mod(angle(X) + 2*pi, 2*pi) / 2 / pi * q), q);
        for aa = 1 : N
        for bb = 1 : N
            L(aa, bb) = A(aa, bb);
        end; end;
    else
        printf('\n * looks like it is not a BH matrix...\n\n', q);
        L = 0;
    end


    if BHF_detector
    %   2025-06-03 BH-fragment detector
        for aa = 1 : N
        r = '';
            for bb = 1 : N
                q = isBH(X(aa, bb), 1024, 'SILENT');
                if q
                    r = strcat(r, 'B');
                else
                    r = strcat(r, '*');
                end
            end; 
            printf("%s\n", r);
        end   
        printf('\n');
    end

end
