function H=sinkhornY8(kMax)
% 2022-08-09
% 2022-08-20 detailed searched for "simple" CHM (small number of HI) with defect=1 or greater than 3
% 2022-12-04
%
% >> sinkhornY8(512 -- 4096)


    addpath ../matrix_tool
    addpath ../CHM
    N = 8;

    while 1
	X = randn(N) + 1j *randn(N);
	Y = sinkhorn(X, kMax);
%        do
%            X = exp(2j*pi*rand(N));
%            Y = sinkhorn(X, kMax);
%            d = ud(Y, "S", 1e-8);
%        until nh(Y) < 1e-13 && n1(Y) < 1e-13 && (d == 0 || d >= 1)


        if nh(Y)<1e-12 && n1(Y)<1e-12
            Y = orderCore(dephase(Y));
            d = ud(Y, "S", 1e-8);
            L = size(getUnique(haagerup(Y), 1e-7),1);
            if L!=10 && L!=70 && L!=74 && L!=130 && L!=170 && L!=1569 && L!=82 && L!=242
                fileName = strcat("Y", int2str(N), "_", int2str(d), "_", int2str(L), "_", datestr(now(), 30), ".dat");
                save(fileName, "Y");
                printf("\nsaved as %s\n", fileName);
            else
                printf("R");
            end
        else
            printf(".");
        end
    end % 1

end



