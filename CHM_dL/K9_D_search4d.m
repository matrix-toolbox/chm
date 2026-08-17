function K9_D_search4d(x0, y0, delta, iiMax)
% 2023-02-16
%
% Searching for inner orbits (defect-4, -16, -12) inside F9_2(zeta).
% Generic defect d(K9) = 2.
% It draws a cloud of points around (x0, y0) with colors representing their defects.
% Parameter delta controls the scope for Gaussian distribution.
%
% >> K9_D_search4d(3, 0, 1, 1024);

    hold all;


    for ii = 1: iiMax
        zeta = x0+randn*delta + 1j*(y0+randn*delta);

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

        %if defect != 2
        %    printf("%d ", defect);
            scatter(real(zeta), imag(zeta), 8, kolor, "filled");
        %end
    end


end


% individual dots with rims:
%
% >> scatter(1,-sqrt(12),64,"y","filled")
% >> scatter(1,-sqrt(12),64,"r")