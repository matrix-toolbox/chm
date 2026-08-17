function H=sinkhornY9(kMax)
% 2022-12-04
%
% >> sinkhornY9(512 -- 4096)


    addpath ../matrix_tool
    addpath ../CHM
    N = 9;
    L0 = 8192;

    while 1
	X = randn(N) + 1j *randn(N);
	Y = sinkhorn(X, kMax);

        if nh(Y)<1e-12 && n1(Y)<1e-12
            isBH(Y, 1000);
            Y = orderCore(dephase(Y));
            d = ud(Y, "S", 1e-8);
            L = size(getUniqueLIMITED(haagerup(Y), L0),1); % filter Y : #L(Y) < 1000
            % printf("#L = %d", L);
	    % #L(H)=105 ==> L9_0 probably...
            if L<L0 && L!=105
                fileName = strcat("Y", int2str(N), "_", int2str(d), "_", int2str(L), "_", datestr(now(), 30), ".dat");
                save(fileName, "*");
                printf("\nsaved as %s\n", fileName);
            else
                printf("R");
            end
        else
            printf(".");
        end
    end % 1

end
