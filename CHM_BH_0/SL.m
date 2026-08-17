function S = SL(M)
% 20170610
% 20200204
% "linear" entropy of matrix M
% entropy is rescaled so that it ranges from 0 to 1
% SL(unitary matrix) = 1
%
% see: _entropy_problem*pdf
%
% the entropy GRM defined is probably proportional to E(U) as E(Ur)=E(S)
% also E(Ut) =E(S) as Utr=Ur S
%
% PRL (2020) E(U) = E_S(U^R)
% "singular entropy" is important for non-unitary matrices (thus we must divide by TrXX' not by the fixed value of dimension
%
%
	%   S = 1 - trace(M*M') / 36/36;
	%S = 1 - sum( (eig(M * M') / size(M, 1)).^2 ); % something's fucky
    
	d = sqrt(size(M,1)); % only square^2 matrices!

	f = M*M'/trace(M*M');
	S = 1 - trace(f*f);
	S = S * (d*d) / (d*d-1);
    
    
    
    
    
    
    


