function [q verdict] = isBH(M, maxExp = 1024, pMode = "VERBOSE")
% ------------------------------------------------------------------------------
% 2017-06-01
%
% Simple check if matrix M is of Butson type.
% 1. real part of entry-wise power should be all ones
% 2. imaginary part ... should be zero
% ------------------------------------------------------------------------------
% Call:
% >> [q verdict] = isBH(S6, 1000, "VERBOSE");
% >> isBH(S6, 1000, "");
% >> isBH(S6, 1000, "VERBOSE");
% ------------------------------------------------------------------------------

    N = size(M, 1);
    verdict = "";

    %b = false;
    q = 0;
    j = 1;
    while j <= maxExp && ~q
        j = j + 1;
        T = M.^j;
        if isequal(real(T), abs(T)) && norm(imag(T), "fro") < 1e-10;
            verdict = sprintf("this is a Butson BH(N=%d, %d)!", N, j);
            %b = true;
            q = j;

            %...................................................................
            % print LOG-form:
            if strcmp(pMode, "VERBOSE")
                LF(M);
            end
            %...................................................................

        end
    end
    if ~q
        verdict = sprintf("probably NOT BH(N=%d, q) for 1 < q < %d...", N, maxExp);
    end
    if strcmp(pMode, "VERBOSE")
        printf("%s\n", verdict);
    end
end

