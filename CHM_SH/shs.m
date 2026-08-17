function shs(N)
% 2023-01-14
% Search for Hermitian (complex) Hadamard matrix using modified Sinkhorn's algorithm.

%   error("BEFORE RUNNING ULOCK PARTICULAR FILTER!");

    while 1
        Y=dephase(sh(exp(2j*pi*rand(N)),4000));
        if nh(Y)<1e-13 && n1(Y)<=1e-13
            A=mod(angle(dephase(Y))+2*pi,2*pi)/2/pi;

            f = 1; % save only matrices with special diagonals: ++++...----
	        print_dA(A);
            for j=1:N/2
                if abs(A(j,j)-0.5)>1e-4, f = 0; end
            end

            L = cL(Y,1e-8);
%           FILTERS:
            if N==8 && L != 18 && L != 54 && L != 102
%	        if N==10 && L != 160 && L != 134 && L != 106 && L != 1002
%	        if N==12 && L != 328 && L != 622 && L != 1902
%	        if f
                d = ud(Y,"S",1e-8);
                printf("N=%d \t d=%d \t L=%d\n", N, d, L);
                save(strcat("SH_", num2str(N), "_", num2str(d), "_L", num2str(L), "_", datestr(now(), 30), ".data"), "A");
            end
        end
    end
end



function print_dA(A)
    N = size(A, 1);
    printf("diagonal = ");
    for j=1:N
        printf("%g ", 2*(2*diag(A)(j)-0.5));
    end
    printf("\n");
end