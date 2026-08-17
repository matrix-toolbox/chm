function Y = DXDsRD(X)
%function Y = DXDsRD(X, p1, p2, p3, p4)

% ------------------------------------------------------------------------------
% 2023-01-02 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% DXD for self R-dual CHM: X = R(X) for X in CHM. Description in "DXDsRG.m".
% ------------------------------------------------------------------------------

    addpath ../matrix_tool
    addpath ../CHM_GITHUB

    r1 = rand;

    w = exp(2j*pi/3);

    N = size(X, 1);
%    while 1
        a1 = rand(1, N);
        a2 = rand(1, N);

        mu = 1;
        Z_OPTIMAL = +Inf;

        while Z_OPTIMAL>7.5e-15
            RESTORE_a1 = a1;
            RESTORE_a2 = a2;
            a1 = mod(a1 + randn(1, N) * mu, 1);
            a2 = mod(a2 + randn(1, N) * mu, 1);
            DL = diag(exp(2j * pi * a1));
            DR = diag(exp(2j * pi * a2));

            DL(1,1) = exp(0j*pi/3);
%            DL(2,2) = exp(0j*pi/3);
%            DL(3,3) = exp(0j*pi/3);
%            DL(4,4) = exp(0j*pi/3);
%            DL(5,5) = exp(0j*pi/3);
%            DL(6,6) = exp(0j*pi/3);
%            DL(7,7) = exp(2j*pi*r1);
%            DL(8,8) = DL(7,7);
%            DL(9,9) = DL(7,7);
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
%            DR(5,5) = exp(2j*pi/9)*exp(2j*pi*p1);
%            DR(6,6) = exp(4j*pi/9)*exp(2j*pi*p2);
%            DR(7,7) = exp(0j*pi/3);
%            DR(8,8) = exp(4j*pi/9)*exp(2j*pi*p3);
%            DR(9,9) = exp(8j*pi/9)*exp(2j*pi*p4);
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
%DR(9,9)/exp(2j*pi*p4)
%diag(DR)

            Y = DL*X*DR;
#SL3(X)
            Z = norm(Y - R(Y),"fro");
            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z * .0031;
                printf("%2.14g\n", Z);
            else
                a1 = RESTORE_a1;
                a2 = RESTORE_a2;
            end
        end




 %   end % while

end

