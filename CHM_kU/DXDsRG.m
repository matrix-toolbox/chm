function Y = DXDsRG(X)
% ------------------------------------------------------------------------------
% 2022-12-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% DXD for self G-dual CHM. Description in "DXD.m".
% ------------------------------------------------------------------------------

    addpath ../matrix_tool
    addpath ../CHM_GITHUB

    y = -1/4+1j*sqrt(15)/4;
    w = exp(2j*pi/3);
    yy = 7/128 + 1j*sqrt((128^2-49)/128^2);

    N = size(X, 1);
    while 1
        a1 = rand(1, N);
        a2 = rand(1, N);

        mu = 1;
        Z_OPTIMAL = +Inf;
#        fileName = strcat("X", num2str(N), "_", datestr(now(), 30), ".data");

        while Z_OPTIMAL>7.5e-15
            RESTORE_a1 = a1;
            RESTORE_a2 = a2;
            a1 = mod(a1 + randn(1, N) * mu, 1);
            a2 = mod(a2 + randn(1, N) * mu, 1);
            DL = diag(exp(2j * pi * a1));
            DR = diag(exp(2j * pi * a2));

%            DL(1,1) = exp(0j*pi/3);
%            DL(2,2) = exp(0j*pi/3);
%            DL(3,3) = exp(0j*pi/3);
%            DL(4,4) = exp(0j*pi/3);
%            DL(5,5) = exp(0j*pi/3);
%            DL(6,6) = exp(0j*pi/3);
%            DL(7,7) = exp(0j*pi/3);
%            DL(8,8) = exp(0j*pi/3);
%            DL(9,9) = exp(0j*pi/3);
%            DL(10,10) = exp(0j*pi/3);
%            DL(11,11) = exp(0j*pi/3);
%            DL(12,12) = exp(0j*pi/3);
%            DL(13,13) = exp(0j*pi/3);
%            DL(14,14) = exp(0j*pi/3);
%            DL(15,15) = exp(0j*pi/3);
%            DL(16,16) = exp(0j*pi/3);
%            DL(9,9)=DL(7,7);
%            DL(7,7) = exp(2j*pi/3);
%            DL(8,8) = exp(4j*pi/3);
%            DL(9,9) = exp(0j*pi/3);
%            DL(10,10) = exp(0j*pi/3);

%            DL(11,11) = exp(0j*pi/3);
%            DL(12,12) = exp(1j*pi/3);
%            DL(13,13) = exp(2j*pi/3);
%            DL(14,14) = exp(3j*pi/3);
%            DL(15,15) = exp(4j*pi/3);
%            DL(16,16) = exp(2j*pi/3);
%diag(DL)
%            DR(1,1) = exp(0j*pi/3);
%            DR(2,2) = exp(0j*pi/3);
%            DR(3,3) = exp(0j*pi/3);
%            DR(4,4) = exp(0j*pi/3);
%            DR(5,5) = exp(4j*pi/3);
%            DR(6,6) = exp(2j*pi/3);
%            DR(7,7) = exp(2j*pi/3);
%            DR(8,8) = exp(4j*pi/3);
%            DR(6,6) = -y;
%            DR(7,7) = -y*y*y;
%            DR(8,8) = yy;
%            DR(9,9) = DR(8,8)*y;
%            DR(7,7) = exp(4j*pi*q2)/exp(2j*pi*q1);
%            DR(8,8) = exp(0j*pi/3);
%            DR(9,9) = exp(0j*pi/3);
%     	    DR(10,10) = exp(2j*pi*q2);
%	        DR(11,11) = exp(3j*pi/3);
%    	    DR(12,12) = 1;
%	        DR(13,13) = -exp(2j*pi*q1)/exp(2j*pi*q2)*exp(1j*pi/3);
%            DR(14,14) = exp(2j*pi*q1)/exp(2j*pi*q2)*exp(1j*pi/3);
%	        DR(15,15) = exp(4j*pi/3);
%     	    DR(16,16) = exp(1j*pi/3);
%DR(8,8)
%diag(DR)

            Y = DL*X*DR;
#SL3(X)
            Z = norm(Y - T(Y),"fro");
            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z * .0031;
                printf("%2.14g\n", Z);
 %               save(fileName, "*");
            else
                a1 = RESTORE_a1;
                a2 = RESTORE_a2;
            end
        end




    end % while

end

