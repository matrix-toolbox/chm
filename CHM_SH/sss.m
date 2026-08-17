function sss(N)
% 2023-01-14
%
% Searching for symmetric (complex) Hadamard matrix using modified Sinkhorn's algorithm.

    while 1
        Y=dephase(ss(exp(2j*pi*rand(N)),4000));
        if nh(Y)<1e-13 && n1(Y)<=1e-14
            L = cL(Y,1e-8);
#            if N==7 && L != 5 && L!= 6 && L != 7
#            if N==8 && L != 10 && L!= 70 && L != 142 && L!=813 && L!=88 && L!=46 && L!=76 && L!=810
#            if N==9 && L != 89 && L!= 201 && L!= 681 && L!=625 && L!=105 && L!=76 && L!=41
#            if N==10 && L != 1071 && L!= 349 && L!= 2071 && L!= 220 && L!=1007 && L!=134 && L!=577 && L!=490
            if N==11 && L > 425 && L < 751 % ... L != 3081 && L!=1457 && L!=1561
                A=mod(angle(dephase(Y))+2*pi,2*pi)/2/pi;
                d = ud(Y,"S",1e-8);
                printf("N=%d \t d=%d \t L=%d\n", N, d, L);
                save(strcat("SS_", num2str(N), "_", num2str(d), "_L", num2str(L), "_", datestr(now(), 30), ".data"), "A");
            end
        end
    end
end

