function P = A36_symmetrize(a0, fi=0)
% -----------------------------------------------------------------------------
% 2023-10-13
% -----------------------------------------------------------------------------
% AME BH(36, 6) can be (easily) symmetrized by appropriate tuning of
% 19 affine parameters.
%
% However, this matrix probably cannot be turn to Hermitian.
% Nor persymmetric...
%
% Input: a0 = vector of non-affine phases.
% It is used only in special cases, otherwise ignored.
% -----------------------------------------------------------------------------
% Symmetric case:
% a = [3 3 3 3 4 1 5 0 3 5 5 4 4 2 4 3 3 3 3];
% or
% a = [2 5 0 0 0 2 1 1 2 2 2 0 0 4 3 2 2 5 0];
% H = A36_BH(a); colormap(prism); imagesc(mod(round(mod(arg(H)+2*pi,2*pi)/2/pi*6),6)); axis square
%
% Unfortunatelly, P is not palindromic nor uniform ... (yet)
%
% 	P=round(LF(A36_BH(a)));
% 	for j=1:36;for k=1:36;if P(j,k)==6,P(j,k)=0;end,end,end
% 	colormap([1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]); imagesc(P); axis square
% -----------------------------------------------------------------------------
% Constant diagonal case:
% It is supposed that for any phase phi in [0, 2pi) one can find such a
% vector P in R^{19} that arg(diag(H(P))) = phi * ones(1, 36).
%
% Example for phi = 0, a 1-parameter family of P's:
% a = rand; P=[5 2 1 a 4 4 2 3 3 0 0 5 5 0 6-a 4 5 2 0];
% H = A36_BH(P); colormap(prism); imagesc(mod(round(mod(arg(H)+2*pi,2*pi)/2/pi*6),6)); axis square
% In particular for:
%    a in {0,1,2,3,4,5}     -->   H in BH(36, 6) with d(B) = 185
%    a in {0,1,2,3,4,5}/k   -->   H in BH(36, 6*k) with d(B) = 185 for k in {1,2,3,4,...}
% -----------------------------------------------------------------------------

	Z_OPTIMAL = +Inf;
	mu = 1;
	P = rand(1, 19);

	while Z_OPTIMAL > 7e-14
		RESTORE_P = P;
		P += randn(1, 19)*mu;
		% --------------------------------------------------------------
		% fixed phases:
		P(1) = 2;
		P(2) = 5;
		P(3) = 0;
		P(4) = 0; % d
		P(5) = 0;
		P(6) = 2;%4;
%%		P(7) = 5;
%%		P(8) = 3;
%%		P(9) = 0; % i
%%		P(10) = 0; % j
%%		P(11) = 0; % k
%%		P(12) = 2; % l
%%		P(13) = 2; % m
%%		P(14) = P(15)-5;
%%%		P(15) = 3*fi/pi;
%%		P(16) = P(15)-1;
%%		P(17) = P(15)+5;
%%		P(18) = P(15)-4;
%%		P(19) = P(15)+3;
		P(19)=0;
		% --------------------------------------------------------------
		X = A36_BH(P);
%		X = A36(a0, P);


		Z = norm(X - transpose(X), "fro"); % <-- symmetric
%		Z = norm(diag(mod(arg(X)+2*pi,2*pi)) - fi*ones(1,36), "fro"); % <-- constant diagonal


		if Z < Z_OPTIMAL
			Z_OPTIMAL = Z
			mu = Z * 0.01;
			save("A36_symmetrize.data", "*");
		else
			P = RESTORE_P;
		end
	end

end

