function Y = HH_12_affine5(x)
% ------------------------------------------------------------------------------
% 2026-09-08 Claude Opus 5
%
% A FIVE-PARAMETER AFFINE FAMILY of complex Hadamard matrices of order 12,
% emanating from the Butson matrix BH(12,4) given by i.^L.
%
%     Y(x) = i.^L .* exp(2i*pi*(x(1)*R1 + ... + x(5)*R5))
%
% Y(x)*Y(x)' = 12*eye(12) identically in x: checked at 200 random points to
% 8.5e-15.  x = 0 gives the Butson matrix itself (defect 15, #Lambda = 4); a
% generic x gives defect 5 and #Lambda = 122.
%
% R2 = R1', and the pair (R1 - R1', R5) spans exactly the HERMITIAN directions,
% a 2-plane.  Restricted to it the family is HH_12_9_32(a,b), of defect 9 and
% #Lambda = 32.  Off that plane the matrices are not Hermitian.
%
% The five directions are the maximal subspace of the 9-dimensional tangent
% space at the Hermitian point on which the moment conditions
%
%     sum_l H_pl conj(H_ql) (R_pl - R_ql)^k = 0,   k = 1..11,   all p < q,
%
% hold; those are exactly equivalent to Y(x) being a CHM for every x, because
% distinct exponentials of the parameter are linearly independent.
% ------------------------------------------------------------------------------

    if nargin < 1, x = zeros(1,5); end
    if numel(x) ~= 5, error('HH_12_affine5: x must have 5 entries'); end

    L = [
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  3  1  2  0  2  1  0  2  2  3;
         0  1  0  0  2  3  1  3  2  2  3  1;
         0  3  0  2  0  2  2  0  2  2  1  0;
         0  2  2  0  0  1  3  3  0  2  2  1;
         0  0  1  2  3  2  1  3  0  0  2  2;
         0  2  3  2  1  3  2  2  0  0  0  1;
         0  3  1  0  1  1  2  2  2  0  3  3;
         0  0  2  2  0  0  0  2  2  2  0  2;
         0  2  2  2  2  0  0  0  2  0  2  0;
         0  2  1  3  2  2  0  1  0  2  0  3;
         0  1  3  0  3  2  3  1  2  0  1  2
    ];

    R1 = [
         0  0  0  0  0  0  0  0  0  0  0  0;
         0 -1 -1 -1  0 -1 -1  0  0  0 -1  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0 -1 -1 -1  0 -1 -1  0  0  0 -1  0;
         0  0  0  0  0  0  0  0  0  0  0  0
    ];

    R2 = [
         0  0  0  0  0  0  0  0  0  0  0  0;
         0 -1  0  0  0  0  0  0  0  0 -1  0;
         0 -1  0  0  0  0  0  0  0  0 -1  0;
         0 -1  0  0  0  0  0  0  0  0 -1  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0 -1  0  0  0  0  0  0  0  0 -1  0;
         0 -1  0  0  0  0  0  0  0  0 -1  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0 -1  0  0  0  0  0  0  0  0 -1  0;
         0  0  0  0  0  0  0  0  0  0  0  0
    ];

    R3 = [
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0 -1 -1 -1  0  0  0 -1;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0 -1  0  0 -1 -1 -1  0  0  0;
         0  0  0 -1  0 -1  0  0 -1  0  0 -1;
         0  0  0 -1  0 -1  0  0 -1  0  0 -1;
         0  0  0  0  0 -1 -1 -1  0  0  0 -1;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0 -1  0  0 -1 -1 -1  0  0  0
    ];

    R4 = [
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0 -1 -1 -1 -1  0  0  0 -1  0 -1  0;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0 -1 -1 -1 -1  0  0  0 -1  0 -1  0;
         0 -1 -1 -1 -1  0  0  0 -1  0 -1  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0 -1 -1 -1 -1  0  0  0 -1  0 -1  0
    ];

    R5 = [
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0 -1 -1  0  0  0  1  0  0  0  1;
         0  1  0  0  0  1  1  1  0  0  1  1;
         0  1  0  0  0  1  1  1  0  0  1  1;
         0  0  0  0  0  1  1  1  0  0  0  1;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0  0 -1 -1 -1  0  0  0 -1  0  0  0;
         0 -1 -1 -1 -1  0  0  0 -1  0 -1  0;
         0  0  0  0  0  1  1  1  0  0  0  1;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0 -1 -1  0  0  0  1  0  0  0  1;
         0 -1 -1 -1 -1  0  0  0 -1  0 -1  0
    ];

    S = x(1)*R1 + x(2)*R2 + x(3)*R3 + x(4)*R4 + x(5)*R5;
    Y = (1i).^L .* exp(2i*pi*S);

end
