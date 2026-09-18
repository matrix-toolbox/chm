function [Y, R] = BH36_6_D_2U_A(P)
% -----------------------------------------------------------------------------
% 2026-08-23
% Wojciech Bruzda | name@cft.edu.pl : name = w.bruzda
% -----------------------------------------------------------------------------
% an affine family of 2-Unitary CHM with 15 parameters stemming from the Dita construction of BH36_6_D_2U
% nine of the parameters act as local diagonal phases, so the family has 6 essential dimensions
% modulo local phases these coincide with the affine family of "A36_BH.m", whose 19 directions split as 13 local phases + the same 6 Dita parameters
% -----------------------------------------------------------------------------
% usage:
%
% >> Y = BH36_6_D_2U_A(P)
%    where P is a vector of ANY phases in R^15 (radians)
%
% >> Y = BH36_6_D_2U_A(rand(1, 15));
%
% in particular
% >> H = BH36_6_D_2U_A(zeros(1, 15)); % is BH36_6_D_2U in BH(36, 6)
% -----------------------------------------------------------------------------
% Y = H o exp(i R) is unimodular and 2-unitary for EVERY real P
% exactly: the two-unitarity conditions are linear in these phases and there are no higher order corrections
% note the convention differs from "A36_BH.m", which uses exp(2i*pi*R/6), i.e., parameters in units of pi / 3
% -----------------------------------------------------------------------------
% origin of R:
%
% in the Dita construction the only deformations not swallowed by a global row or column rescaling are the column scalings inside the twelve order-three blocks
%
% F --> F * diag(exp(i psi(J, 0)), exp(i psi(J, 1)), exp(i psi(J, 2)))
%
% one diagonal per block
% [ulled back to the four-leg indexing of BH36_6_D_2U with rows 6 * a + b + 1  and columns 6 * c + d + 1
% these are the phase arrays
%
% R(6 * a + b + 1, 6 * c + d + 1) = psi(J, be)
% J = 2 * a + mod(b, 2)
% be = mod(d, 3)
%
% a 36-dimensional space
%
% fourteen of those dimensions are plain row/column rescalings, which leaves Dita's (12 - 1) * (3 - 1) = 22 essential parameters
%
% imposing 2-unitarity on the 22 is an exactly LINEAR condition and retains 15 of the 36 phases
% the surviving solutions are exactly those with
%
% 1. psi(a, 0, be) - psi(a, 1, be) = D(be + s(a)) for ONE function D on Z_3
%    with s(a) = mod(a + 1, 3)
%    the same s(a) as the row rescaling of BH36_6_D_2U
% 2. psi(a, 1, be) - psi(a + 3, 1, be) independent of be, for a = 0, 1, 2
%
% where psi(a, b2, be) denotes psi(2 * a + b2, be)
%
% Free data: c (3 x 3), k (3), D (3)
% hence 9 + 3 + 3 = 15
% explicitly the 12 x 3 table psi(J, be) reads
%
%                 be =     0              1              2
% J =  0  (a = 0, b2 = 0)  c11 + D2       c12 + D3       c13 + D1
% J =  1  (a = 0, b2 = 1)  c11            c12            c13
% J =  2  (a = 1, b2 = 0)  c21 + D3       c22 + D1       c23 + D2
% J =  3  (a = 1, b2 = 1)  c21            c22            c23
% J =  4  (a = 2, b2 = 0)  c31 + D1       c32 + D2       c33 + D3
% J =  5  (a = 2, b2 = 1)  c31            c32            c33
% J =  6  (a = 3, b2 = 0)  c11 + k1 + D2  c12 + k1 + D3  c13 + k1 + D1
% J =  7  (a = 3, b2 = 1)  c11 + k1       c12 + k1       c13 + k1
% J =  8  (a = 4, b2 = 0)  c21 + k2 + D3  c22 + k2 + D1  c23 + k2 + D2
% J =  9  (a = 4, b2 = 1)  c21 + k2       c22 + k2       c23 + k2
% J = 10  (a = 5, b2 = 0)  c31 + k3 + D1  c32 + k3 + D2  c33 + k3 + D3
% J = 11  (a = 5, b2 = 1)  c31 + k3       c32 + k3       c33 + k3
%
% and R below is that table lifted to 36 x 36
%
% R is constant in c and 3-periodic in d, so every row is one three-entry pattern repeated 12 times
% rows with the same a and the same parity of b coincide
% -----------------------------------------------------------------------------

	if nargin < 1 || isempty(P)
		P = zeros(1, 15);
	end
	if numel(P) != 15 || !isreal(P)
		error("parameter P must be a real vector of 15 phases!");
	end

	c11 = P(1);
	c12 = P(2);
	c13 = P(3);
	c21 = P(4);
	c22 = P(5);
	c23 = P(6);
	c31 = P(7);
	c32 = P(8);
	c33 = P(9);
	k1  = P(10);
	k2  = P(11);
	k3  = P(12);
	D1  = P(13);
	D2  = P(14);
	D3  = P(15);

	R = [
		c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1;
		c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13;
		c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1;
		c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13;
		c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1      c11+D2    c12+D3    c13+D1    c11+D2    c12+D3    c13+D1;
		c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13         c11       c12       c13       c11       c12       c13;

		c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2;
		c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23;
		c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2;
		c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23;
		c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2      c21+D3    c22+D1    c23+D2    c21+D3    c22+D1    c23+D2;
		c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23         c21       c22       c23       c21       c22       c23;

		c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3;
		c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33;
		c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3;
		c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33;
		c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3      c31+D1    c32+D2    c33+D3    c31+D1    c32+D2    c33+D3;
		c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33         c31       c32       c33       c31       c32       c33;

		c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1;
		c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1;
		c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1;
		c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1;
		c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1   c11+k1+D2 c12+k1+D3 c13+k1+D1 c11+k1+D2 c12+k1+D3 c13+k1+D1;
		c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1      c11+k1    c12+k1    c13+k1    c11+k1    c12+k1    c13+k1;

		c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2;
		c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2;
		c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2;
		c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2;
		c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2   c21+k2+D3 c22+k2+D1 c23+k2+D2 c21+k2+D3 c22+k2+D1 c23+k2+D2;
		c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2      c21+k2    c22+k2    c23+k2    c21+k2    c22+k2    c23+k2;

		c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3;
		c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3;
		c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3;
		c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3;
		c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3   c31+k3+D1 c32+k3+D2 c33+k3+D3 c31+k3+D1 c32+k3+D2 c33+k3+D3;
		c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3      c31+k3    c32+k3    c33+k3    c31+k3    c32+k3    c33+k3;
	];

	Y = BH36_6_D_2U() .* exp(1j * R);

end
