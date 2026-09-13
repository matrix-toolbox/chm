function Lambda = haagerup(M, d_mode="SILENT")
% ------------------------------------------------------------------------------
% 2016-06-04
% 2017-12-12
% Haagerup LAMBDA invariants set
%
% 2021-09-19
% ------------------------------------------------------------------------------

    if (!strcmp(d_mode, "SILENT"))
        disp(sprintf('* wait: calculating Haagerup invariants...'));
    end

    N = size(M, 1); % used only for square (CHM) matrices
    % Lambda  = zeros(N^4, 1);
    % mu = 1;
    Lambda = [];
    for a = 1 : N
        if (!strcmp(d_mode, "SILENT"))
        printf('.')
    end
    for b = 1 : N
    for c = 1 : N
    for d = 1 : N
    %   Lambda(mu) = M(a, b) * M(c, d) * conj(M(a, d)) * conj(M(c, b));
    %   mu = mu + 1;
        Lambda = [Lambda; M(a, b) * M(c, d) * conj(M(a, d)) * conj(M(c, b))];
    end; end; end; end

    if (!strcmp(d_mode, "SILENT"))
        printf('\n');
    end

end

