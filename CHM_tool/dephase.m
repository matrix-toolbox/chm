function K = dephase(H)
% 20200831
% dephased (normalized) form of H with 1st row and 1st column of ones
%
% 1 1 1 1 1 1 ... 1
% 1 . . . . . ... .
% 1 . . . . . ... .
% 1 CORE(d-1)x(d-1)
% 1 . . . . . ... .
% 1 . . . . . ... .
% : : : : : :     :
% 1 . . . . . ... .

    H = H * diag(1 ./ H(1,:));
    K = diag(1 ./ H(:, 1)) * H;

    K(1,:) = ones(1, size(H, 1)); % additional cleaning to get rid of annoying zeros: 0.9517898237498213e-16
    K(:,1) = ones(size(H, 1), 1);

end
