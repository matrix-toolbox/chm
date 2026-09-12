function Y = DS_16_9_L(u, v)
% ------------------------------------------------------------------------------
% 2026-09-05 Claude Opus 5 Max ==> (semi-)analytic structure from numerical pattern
%
% nine-parameter AFFINE family of dephased complex Hadamard matrices of order 16
% containing the numerical point DS_16_9_L_20240117T094601.data
%
%   Y = DS_16_9_L()        % the stored numerical point
%   Y = DS_16_9_L(w)       % w    1x5 : symmetric slice, Y == Y.'  (5 parameters)
%   Y = DS_16_9_L(u, v)    % u,v  1x5 : the full family            (9 parameters)
%
% STRUCTURE
%   Rows/columns 1,2,3,12,15,16 are frozen.  The remaining ten split into five
%   PAIRS, which are the orbits of the exact involutive automorphism
%       pi = (4 9)(5 8)(6 13)(7 14)(10 11)      % A(pi(i),pi(j)) == A(i,j)
%
%       P1 = {4,9}   P2 = {5,8}   P3 = {6,13}   P4 = {7,14}   P5 = {10,11}
%
%   These cut the active part into a 5x5 grid of 2x2 cells.  The five diagonal
%   cells are frozen; each off-diagonal cell (Pi,Pj) carries ONE free phase
%
%       x(i,j) = u(i) + v(j)          (mod 1),      i ~= j
%
%   so the family has 5 + 5 - 1 = 9 real parameters ((u,v) -> (u+t, v-t) acts
%   trivially).  9 equals the defect of every generic member, hence this family
%   is MAXIMAL: nothing is missing.
%
%   The base point (u = v = 0) is a symmetric BH(16,6) - all phases are
%   multiples of 1/6.  It is more degenerate than the family (defect 13).
%
% LETTER FORM (the a,b,c,d,e,f,g,h,k,m of the earlier attempt)
%   Each letter is one unordered pair of the five blocks:
%       a={1,2} b={1,3} c={1,4} f={1,5} d={2,3} e={2,4} g={2,5} k={3,4} h={3,5} m={4,5}
%   On the symmetric slice (u = v = w) they are  letter{i,j} = base + w(i) + w(j),
%   so ONLY FIVE of the ten are free.  Taking a,b,c,d,f as free:
%
%       e = c - b + d + 1/3
%       g = f - b + d + 1/6
%       k = c + d - a + 2/3
%       h = f + d - a + 1/3
%       m = c + f + d - a - b
%
%   (the first two are the relations found earlier; the last three were missing,
%   which is why letting a,b,d,e,g,h,k,m vary independently broke orthogonality)
% ------------------------------------------------------------------------------
 
%   exponent matrix of the BH(16,6) base point: phases are E/6
    E = [
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 3 2 5 1 5 1 1 5 3 3 2 5 1 4 4;
        0 2 3 3 1 1 5 1 3 5 5 4 1 5 2 4;
        0 5 3 3 3 0 0 0 3 3 0 1 3 3 5 1;
        0 1 1 3 3 3 5 3 0 4 1 5 0 2 5 3;
        0 5 1 0 3 3 4 0 3 5 2 3 3 1 1 5;
        0 1 5 0 5 4 3 2 3 3 0 1 1 3 3 5;
        0 1 1 0 3 0 2 3 3 1 4 5 3 5 5 3;
        0 5 3 3 0 3 3 3 3 0 3 1 0 0 5 1;
        0 3 5 3 4 5 3 1 0 3 3 5 2 0 1 1;
        0 3 5 0 1 2 0 4 3 3 3 5 5 3 1 1;
        0 2 4 1 5 3 1 5 1 5 5 3 3 1 4 2;
        0 5 1 3 0 3 1 3 0 2 5 3 3 4 1 5;
        0 1 5 3 2 1 3 5 0 0 3 1 4 3 3 5;
        0 4 2 5 5 1 3 5 5 1 1 4 1 3 3 2;
        0 4 4 1 3 5 5 3 1 1 1 2 5 5 2 3
    ];
 
    PAIRS = [4 9; 5 8; 6 13; 7 14; 10 11];
 
    if nargin == 0
%       the stored numerical point (symmetric slice)
        u = [0.90277718490137937 0.27034682639618701 0.09823773963693141 0.45342732559660193 0.57478637156408263];
        v = u;
    elseif nargin == 1
        v = u;
    end
 
    if numel(u) ~= 5 || numel(v) ~= 5
        error('DS_16_9_L: u and v must each have 5 entries.');
    end
    if ~isnumeric(u) || ~isnumeric(v) || any(~isreal([u(:); v(:)])) ...
                                      || any(~isfinite([u(:); v(:)]))
        error('DS_16_9_L: parameters must be real and finite.');
    end
 
    R = zeros(16, 16);
    for i = 1:5
        for j = 1:5
            if i == j
                continue;
            end
            R(PAIRS(i,:), PAIRS(j,:)) = u(i) + v(j);
        end
    end
 
    Y = exp(2j * pi * mod(E / 6 + R, 1));
end
