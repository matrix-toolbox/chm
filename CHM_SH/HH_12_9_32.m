function Y = HH_12_9_32(a, b)
% ------------------------------------------------------------------------------
% 2026-09-08 Claude Opus 5
%
% Hermitian CHM(12) with d = 9 and L = 32
%
% This is not an isolated matrix: it is a point of a TWO-PARAMETER AFFINE FAMILY
% of Hermitian complex Hadamard matrices.  The individual entries of any single
% member are therefore not algebraic numbers -- the family is the analytic
% object, and (a,b) are free real moduli.
%
%     Y(a,b) = i.^L .* exp(2i*pi*(a*P + b*Q))
%
% Y(a,b)*Y(a,b)' = 12*eye(12) holds identically in a and b (verified to 60
% digits at irrational values), and Y = Y' for every a,b because P and Q are
% antisymmetric.
%
%   a = b = 0   ->  i.^L, an exact BH(12,4): every entry a 4th root of unity,
%                   Hermitian, defect 15, #Lambda = 4.  The family is the
%                   affine deformation of that Butson matrix.
%   generic a,b ->  defect 9, #Lambda = 32.
%
% Different (a,b) give monomially INEQUIVALENT matrices, so the parameters are
% genuine moduli and not a re-gauging.
%
% P is a pure off-diagonal block: P = 1_(TxS) - 1_(SxT) with S = {2,11} and
% T = {3,4,6,7} in 1-based indexing.
%
% The defaults reproduce the numerically found HH_12_9_32_20260908T103306.m to
% 7.3e-13.  This Hermitian plane sits inside a 5-parameter affine family of
% general CHM(12) -- see HH_12_affine5.m.
% ------------------------------------------------------------------------------

    if nargin < 1, a = 0.050655317368224415; end
    if nargin < 2, b = 0.078661517462416272; end

    % quarter-turn exponents: the BH(12,4) skeleton
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

    % first antisymmetric affine direction (a pure off-diagonal block)
    P = [
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0 -1 -1  0 -1 -1  0  0  0  0  0;
         0  1  0  0  0  0  0  0  0  0  1  0;
         0  1  0  0  0  0  0  0  0  0  1  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  1  0  0  0  0  0  0  0  0  1  0;
         0  1  0  0  0  0  0  0  0  0  1  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0;
         0  0 -1 -1  0 -1 -1  0  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  0
    ];

    % second antisymmetric affine direction
    Q = [
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

    Y = (1i).^L .* exp(2i*pi*(a*P + b*Q));

end
