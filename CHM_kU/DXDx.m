function Y = DXDx(X)
% ------------------------------------------------------------------------------
% 2022-12-21 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2022-12-26
%
% Check whether Y = D1*X*D2 can be a tensor product for two diagonals Dj
% and X being some CHM of order N=9 or 16.
% ------------------------------------------------------------------------------

#p1 = rand;
#p2 = rand;
    addpath ../matrix_tool
    addpath ../CHM_GITHUB

    N = size(X, 1);
    d = floor(sqrt(N));
    while 1
        a1 = rand(1, N);
        a2 = rand(1, N);
        mu = 1;
        Z_OPTIMAL = +Inf;
%        fileName = strcat("X", num2str(N), "_", datestr(now(), 30), ".data");


        while Z_OPTIMAL>7.5e-15
            RESTORE_a1 = a1;
            RESTORE_a2 = a2;
            a1 = mod(a1 + randn(1, N) * mu, 1);
            a2 = mod(a2 + randn(1, N) * mu, 1);
            DL = diag(exp(2j * pi * a1));
            DR = diag(exp(2j * pi * a2));

%            DL(1,1) = 1;
%            DL(2,2) = 1;
%            DL(3,3) = 1;
%            DL(4,4) = 1;
%            DL(5,5) = 1;
%            DL(6,6) = 1;
%            DL(7,7) = exp(2j*pi*p1);
%            DL(8,8) = DL(7,7);
%            DL(9,9) = exp(2j*pi/3);
%            DL(10,10) = exp(2j*pi/3);
%            DL(11,11) = -1;
%            DL(12,12) = -1;
%            DL(13,13) = -1;
%            DL(14,14) = 1;
%            DL(15,15) = -1;
%            DL(16,16) = 1;
%aaa=DL(10,10)/exp(2j*pi/3)
%Z_OPTIMAL
%ddl=diag(DL)
%            DR(1,1) = 1;
%            DR(2,2) = 1;
%            DR(3,3) = 1;
%            DR(4,4) = 1;
%            DR(5,5) = 1;
%            DR(6,6) = 1j;
%            DR(7,7) = 1;
%            DR(8,8) = 1j;
%            DR(9,9) = exp(2j*pi*p2);
%            DR(10,10) = exp(2j*pi*p2);
%            DR(11,11) = exp(2j*pi*p2);
%            DR(12,12) = exp(2j*pi*p2);
%            DR(13,13) = 1;
%            DR(14,14) = 1j;
%            DR(15,15) = 1;
%            DR(16,16) = 1j;
%aa1=DR(7,7)*exp(1j*pi*p3)/exp(1j*pi*p2)
%aa2=exp(-2j*pi*p2)
%aa2=DR(15,15)
%S_OPTIMAL
%ddr=diag(DR)

            Y = DL*X*DR;
            Z = isTensor(Y,d,d,d,d); %abs(SL(Y) - 1) + abs(SL(R(Y)) - 1) + abs(SL(T(Y)) - 1);
            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z * .000021;
                printf("%2.14g\n", Z);
%               save(fileName, "*");
            else
                a1 = RESTORE_a1;
                a2 = RESTORE_a2;
            end
        end




    end % while

end

