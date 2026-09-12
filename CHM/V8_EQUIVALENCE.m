function V8_EQUIVALENCE
% ------------------------------------------------------------------------------
% 2021-09-19 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Proof that original matrix of Veit V_8^{(0)} % is equivalent
% to the one obtained from the sequence of V8 of isolated non-Butson CHM.
%
% NOTE THAT THIS IS ONLY ONE OF MANY POSSIBLE CONFIGURATIONS
% FOR a, b and c (\mp, c.c.).
%
% Similarly, one proves that A_8^{(0)} is equivalent to another
% solution of V8 (sequence).
% -----------------------------------------------------------------------------

	clc;
	% analytic values for a, b and c can be taken from V8_ANALYTIC.[m|nb]

	a =  0.650891015304204460 + 0.75917118372358006*1j;
	b =  0.194170403372279250 - 0.98096781519795351*1j;
	c = -0.032552288838649801 - 0.99947003381360311*1j;

        VEIT_ORIGINAL_UNDEPHASED=[
        	-1  -1   b   b   c   c   a   a;
        	-1   b  -1   c   b   a   c  -a;
        	 b  -1   c  -1   a   b  -a   c;
        	 b   c  -1   a  -1  -a   b  -c;
        	 c   b   a  -1  -a  -1  -c   b;
        	 c   a   b  -a  -1  -c  -1  -b;
        	 a   c  -a   b  -c  -1  -b  -1;
        	 a  -a   c  -c   b  -b  -1   1;
        ];

	LD = diag([-1, -1, 1/b, 1/b, 1/c, 1/c, 1/a, 1/a]);
	RD = diag([1, 1, -1/b, -1/b, -1/c, -1/c, -1/a, -1/a]);
	VEIT_ORIGINAL_DEPHASED = LD * VEIT_ORIGINAL_UNDEPHASED * RD

	load V8_V8_0.dat % get the matrix from the sequence V_8
	V8_SEQUENCE = dephase(H);
	V8_FROM_SEQUENCE_BEFORE_PERMUTED = V8_SEQUENCE
	
	LP = [ % permute rows
		1 0 0 0 0 0 0 0;
		0 0 0 0 0 0 0 1;
		0 1 0 0 0 0 0 0;
		0 0 0 0 0 0 1 0;
		0 0 1 0 0 0 0 0;
		0 0 0 0 0 1 0 0;
		0 0 0 1 0 0 0 0;
		0 0 0 0 1 0 0 0;
	];
	RP = [ % permute columns
		1 0 0 0 0 0 0 0;
		0 1 0 0 0 0 0 0;
		0 0 0 1 0 0 0 0;
		0 0 0 0 0 1 0 0;
		0 0 0 0 0 0 0 1;
		0 0 0 0 0 0 1 0;
		0 0 0 0 1 0 0 0;
		0 0 1 0 0 0 0 0;
	];

	V8_FROM_SEQUENCE_AFTER_PERMUTED = LP*V8_SEQUENCE*RP
		
	DIFFERENCE = round(VEIT_ORIGINAL_DEPHASED - V8_FROM_SEQUENCE_AFTER_PERMUTED)

end
