function name = bc_name(prefix, H, A)
% ------------------------------------------------------------------------------
% 2026-09-13 Wojciech Bruzda, name[at]uj.edu.pl : name = w.bruzda, https://matrix-toolbox.github.io/chm
% ------------------------------------------------------------------------------
% File name for a block-circulant solution, following the convention used
% throughout the Catalog:
%
%     <prefix>_<N>_<defect>_<#Lambda>.data       e.g.  VH_10_0_109.data
%
% Called by the M-files that BC_L.sh and C_V.sh generate.  H is the matrix,
% A its phase form; only A is stored.  Both invariants are expensive, so they
% are computed only while that is affordable and the name degrades gracefully:
%
%     N <=  40   <prefix>_<N>_<d>_<L>.data
%     N <=  60   <prefix>_<N>_<d>.data           (#Lambda too slow)
%     N >   60   <prefix>_<N>.data               (both too slow)
%
% A letter is appended when the name is taken, so a second solution with the
% same invariants does not overwrite the first: VH_14_0_687A.data, ...B, ...
%
% >> save(bc_name('VH', V, A), 'A');
% ------------------------------------------------------------------------------

    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'CHM_BH_0'));
    N = size(H, 1);

    name = sprintf('%s_%d', prefix, N);
    if N <= 60
        name = sprintf('%s_%d', name, ud(H, 'S', 1e-8));
        if N <= 40
            name = sprintf('%s_%d', name, nLambda(H));
        end
    end

    stem = name;
    name = [stem '.data'];
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    k = 0;
    % dir(), not exist(): exist() also searches the load path, so a copy of the
    % Catalog on the path would make every run look like a collision
    while ~isempty(dir(name))
        k++;
        if k > numel(letters)
            error('bc_name: too many solutions named %s', stem);
        end
        name = sprintf('%s%s.data', stem, letters(k));
    end

end


function n = nLambda(H)
% Cardinality of the Haagerup set, #Lambda(H) = |{ H_ij H_kl conj(H_il) conj(H_kj) }|
    N = size(H, 1);
    L = zeros(N^4, 1);
    t = 0;
    for i = 1:N, for j = 1:N, for k = 1:N, for l = 1:N
        t++;
        L(t) = H(i, j) * H(k, l) * conj(H(i, l)) * conj(H(k, j));
    end, end, end, end
    n = rows(unique(round([real(L) imag(L)] * 1e6) / 1e6, 'rows'));
end
