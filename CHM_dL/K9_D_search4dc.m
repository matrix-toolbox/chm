function K9_D_search4dc(delta, iiMax, orbit)
% 2023-02-16
%
% Draw a random point from a given (sub)orbit (see: K9_suborbit.m)
% and check its vicinity for the defect.
%
% >> k = 1; % 2 3 4
% >> K9_D_search4dc(0.1, 10, 1);

    hold all;


    for ii = 1: 4096

        t = rand;
        if orbit == 1
            zeta = 1 + 2*exp(2j*pi*t); % CR
        elseif orbit == 2
            zeta = -1 + 2*exp(2j*pi*t); % CL
        elseif orbit == 3
            zeta = 1j*sqrt(3) + 2*exp(2j*pi*t); % CU
        elseif orbit == 4
            zeta = -1j*sqrt(3) + 2*exp(2j*pi*t); % CD
        end
        
        K9_search4d(real(zeta), imag(zeta), delta, iiMax);

    end



end
