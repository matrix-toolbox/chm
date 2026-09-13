function C = core(H)
% 2025-06-23
% get the core of H

    N = size(H, 1);
    C = H(2:N, 2:N);

end
