function Y = DXD(X, p1, p2, p3, p4)
%function Y = DXD(X)

% ------------------------------------------------------------------------------
% 2022-12-21 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2022-12-26
%
% Check whether Y = D1*X*D2 can be 2-Unitary for two diagonals Dj and X being
% some CHM of order N = 9 or 16.
% Two (unitary) diagonal matrices have arbitrary entries.
% Random walk procedure tries to optimize them so that Y is 2U-CHM.
% When particular entries of DL or DR are fixed, one can recover analytical
% structure of Y.
% ------------------------------------------------------------------------------
% Sometimes entries of DL (DR) depend on parameters pj from X(pj).
% To write them correctly, either use Mathematica to solve unitary constraints
% for R(X) and T(X) (bad idea...), or... just guess them!
% In many cases they are simple functions of pj like 1/pj, w*pj^2, pj^2/pk, etc.
%
% When there are two or more parameters pj, it is a good idea to set to zero
% all but one to check which entries of Dx depend on it. Checking in this way
% all parameters, one easily discovers which Dx(k,k) is a function of
% (a) particular parameter(s) pj.
%
% Or, even better, set each single pj to exp(2j*pi/7) and check the output.
% For example:
%
% >> p1=1/7; p2=0;   p3=0;   DXD(Y16_d12_L90(p1, p2, p3), p1, p2, p3)
% >> p1=0;   p2=1/7; p3=0;   DXD(Y16_d12_L90(p1, p2, p3), p1, p2, p3)
% >> p1=0;   p2=0;   p3=1/7; DXD(Y16_d12_L90(p1, p2, p3), p1, p2, p3)
%
% Knowing numerical value of exp(2j*pi/7) or one another "exotic" root
% of unity, one immediately recognizes (all) internal dependencies!
% However, this works only for affine families!
% For non-affine one (if there is any?) this might be much more complicated,
% and probably analytically intractable...
% ------------------------------------------------------------------------------
%
% Usage:
% >> DXD(Z9(rand, rand))
% >> DXD(kron(F3, F3))
% >> DXD(B9q3(1))
%
% Works for many 9- and 16-dimensional CHM which are not isolated.
% In particular:
%
% Y = DL * kron(F2, F4(0), F2) * DR for some DL, DR
% Y = DL * kron(F5, F5) * DL
% ******************************************************************************
% DOES NOT work (so far) for:
% N=9
%     B9
%     N9_0
%     S9_0 (smallest Z)
%     Y9_0 ("block circulant")
%     some BH(9, q), in particular BH(9, 10)
%
% N=16
%     kron(F2, F2, F2, F2)
%     kron(F2, S8_4(0,0,0,0))
%     kron(S8_4(0,0,0,0), F2)
%     L16_0 <-- this supports the conjecture that
%               isolated CHM do not provide 2UCHM
%     F16_17(a) : a = rand(1,17) || zeros(1,17) || ones(1,17)/3
%     S16_11(a) : a = ...
%     H16x for x = A, B, C, D, E
%     no matrix taken at random from BH-16-4 (from BUTSON Home!)
%
% N=36
%     kron(F6_2(1/2,1/3), F6_2(1/3,1/2)) ---> 0.0131
% ******************************************************************************

r1 = rand;
r2 = rand;
    addpath ../matrix_tool
    addpath ../CHM_GITHUB

    N = size(X, 1);
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

#            DL(1,1) = 1;
#            DL(2,2) = 1;
#            DL(3,3) = 1;
#            DL(4,4) = 1;
#            DL(5,5) = exp(4j*pi/3);
#            DL(6,6) = exp(2j*pi/3);
#            DL(7,7) = -1j/exp(2j*pi/3);
#            DL(8,8) = -1j;
#            DL(9,9) = -1j*exp(2j*pi/3);

%            DL(10,10) = exp(2j*pi/3);
%            DL(11,11) = -1;
%            DL(12,12) = -1;
%            DL(13,13) = -1;
%            DL(14,14) = 1;
%            DL(15,15) = -1;
%            DL(16,16) = 1;
%aaa=DL(9,9)
%S_OPTIMAL
%ddl=diag(DL)
#            DR(1,1) = 1;
#            DR(2,2) = 1;
#            DR(3,3) = 1;
#            DR(4,4) = 1;
#            DR(5,5) = exp(2j*pi/3);
#            DR(6,6) = exp(4j*pi/3);
#            DR(7,7) = exp(2j*pi*r1);
#            DR(8,8) = exp(2j*pi*r1)/exp(2j*pi/3);
#            DR(9,9) = exp(2j*pi*r1)/exp(4j*pi/3);

%            DR(10,10) = exp(2j*pi*p2);
%            DR(11,11) = exp(2j*pi*p2);
%            DR(12,12) = exp(2j*pi*p2);
%            DR(13,13) = 1;
%            DR(14,14) = 1j;
%            DR(15,15) = 1;
%            DR(16,16) = 1j;
%aa1=DR(5,5)
%aa2=exp(-2j*pi*p2)
%aa2=DR(15,15)
%S_OPTIMAL
%ddr=diag(DR)
            Y = DL*X*DR;
#SL3(Y)
            Z = abs(SL(Y) - 1) + abs(SL(R(Y)) - 1) + abs(SL(T(Y)) - 1);
            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z * 2.31;
                printf("%2.14g\n", Z);
%               save(fileName, "*");
            else
                a1 = RESTORE_a1;
                a2 = RESTORE_a2;
            end
        end
    end % while

end

