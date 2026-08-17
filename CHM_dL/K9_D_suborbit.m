function K9_D_suborbit(orbit)
% 2023-02-16
%
% Draw one of four orbits == circles inside D.
%
% >> k = 1; % 2 3 4
% >> K9_D_suborbit(k);
%

    hold all;



    for t = 0:0.001:1
        if orbit == 1
            zeta = 1 + 2*exp(2j*pi*t); % CR
        elseif orbit == 2
            zeta = -1 + 2*exp(2j*pi*t); % CL
        elseif orbit == 3
            zeta = 1j*sqrt(3) + 2*exp(2j*pi*t); % CU
        elseif orbit == 4
            zeta = -1j*sqrt(3) + 2*exp(2j*pi*t); % CD
        end
        

        H = K9_2z(zeta);
        if H == -1, continue, end

        defect = ud(H, "S", 1e-8);
        if defect == 2
            kolor = "r";
        elseif defect == 4
            kolor = "g";
        elseif defect == 10
            kolor = "b";
        elseif defect == 12
            kolor = "c";
        elseif defect == 16
            kolor = "m";
        end
        scatter(real(zeta), imag(zeta), 4, kolor, "filled")
    end



end


