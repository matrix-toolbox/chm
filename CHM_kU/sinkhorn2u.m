function Y = sinkhorn2u(N, kMax)
% ------------------------------------------------------------------------------
% 2022-11-15 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Modified Sinkhorn algorithm to search for 2-unitary Hadamards.
%
% Usage:
%
% >> sinkhorn2u(9,2^10)
% >> sinkhorn2u(16,2^13)
% >> sinkhorn2u(25,2^16)
%
% For N=9 or N=16 all numerically found matrices with low #Lambda can be
% transformed into affine families...
%
%        procedure:
%        1. sinkhorn2u(N, kMax) ---> YN_d.._L.._TIMESTAMP.data
%        2. DXD.m to recover DL and DR (usually internal structure of Y should not be changed!)
%        3. ---> YN_dx_Lx(pj).m which is 2U: SL(Y..., "minimal") = (1, 1, 1)
%
%        "sinkhorn" outputs Y : SL3(Y, "minimal") = (1, 1, 1)
%
%        for N = 16 sinkhorn2u returns many matrices with various of defects.
%        and very different #Lambdas -- this might be numerical discrepancies!
%        All of these matrices are recovered by analytical affine families: Y9_dx_Lx(pj).m
% ------------------------------------------------------------------------------

    addpath ../AME46
    addpath ../CHM
    addpath ../matrix_tool

    while 1
        % X = randn(N) + 1j * randn(N);
        G = randn(N) + 1j*randn(N);
        %X = P36_tilde * expm(1j * 1e-8 * (G+G')/2);
        X = PM(N) * expm(1j * 10^-randi(10) * (G+G')/2);

        k = 0;
        while (nh(X)>=1e-14 || n1(X)>=5e-16) && (k<=kMax)
            X = X ./ abs(X);
            X = R(T(PD(X)));
            k++;
        end

        Y = X*sqrt(N); % |H| = 1, H*H' = N*I

        S1 = SL(Y);
        S2 = SL(R(Y));
        S3 = SL(T(Y));
        norm1 = nh(Y);
        norm2 = n1(Y);

	KNOWN = [ % { N, d, #L }
	    9, 10, 21;
	    9, 16, 3;
	    9, 2, 89;
	    9, 4, 21;
	    9, 4, 27;
	    9, 4, 3;
            9, 4, 39;
	    9, 4, 57;
	    16, 21, 370;
	    16, 19, 90;
	    16, 17, 306;
	    16, 17, 918;
	];

	if abs(S1-1)<1e-13 && abs(S2-1)<1e-13 && abs(S3-1)<1e-13 && norm1<1e-13 && norm2<1e-13
		d = ud(Y, "S", 1e-8);
%%		L = size(getUnique(haagerup(Y), 1e-7), 1);

%%		f = 0;
%%		for j=1:size(KNOWN, 1)
%%			k = KNOWN(j, :);
%%			if (k(1)==N) && (k(2)==d) && (k(3)==L)
%%				f = 1;
%%				break;
%%			end
%%		end

%%		if !f
%%			fileName = strcat("Y", int2str(N), "_d", int2str(d), "_L", int2str(L), "_", datestr(now(), 30), ".dat");
			fileName = strcat("Y", int2str(N), "_d", int2str(d), "_L_", datestr(now(), 30), ".data");
			A = mod(angle(Y)+2*pi,2*pi)/2/pi;
			save(fileName, "A");
			printf("\n SL3 = (%g %g %g) \t nh = %g \t n1 = %g ", S1, S2, S3, norm1, norm2);
			printf("\t ----> saved as %s", fileName);
%%		else
%%			printf("."); % known triplet
%%		end
	else
		printf("x"); % fail -- algorithm does not converge
	end

    end

end

