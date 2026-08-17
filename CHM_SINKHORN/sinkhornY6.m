function Y = sinkhornY6(kMax)
% 2022-08-09
% 2022-08-20
%
% Detailed searched for a "simple" CHM (small number of HI) with defect = 4 FAILS.
% It only finds Y : #L(Y) = 3 or 451. The latter one without any simple internal dependencies.
% Sometimes a symmetric matrix is observed, but only for Y : #L(Y) = 3 (spectral one).

    N = 6;
    while(1)
        do
            X = exp(2j*pi*randn(N));
            Y = sinkhorn(X, kMax);
            defect = ud(Y,"S",1e-8);
        until nh(Y) < 1e-14 && n1(Y) < 1e-14 && defect > 0 % && (defect == 0 || defect >= 1)

        Y = dephase(Y);
        cSH(Y);
        cc = sfc(Y, "SILENT");
        if cc < 25
            fileName = strcat('Y', int2str(N), '_', int2str(defect), '_L', int2str(cL(Y,1e-8)), '_', datestr(now(), 30), '.dat');
            printf("save: %s\n", fileName);
            A = mod(angle(Y)+2*pi, 2*pi)/2/pi;
            save(fileName, "A");
        end
    end
end
