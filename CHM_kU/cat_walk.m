function U = cat_walk(d, a0, b0, c0, mu0)
% ------------------------------------------------------------------------------
% 2022-06-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2022-11-21
%
% >> cat_walk(d, a0, b0, c0, mu0) # N = d^2
%
% output is usually very inaccurate, eg.:
%    a=-2.0119        b=-7.01925      c=-4,   Z=0.000793754   [#1]
%
% after "fine tuning" one recovers 2-Unitary CHM:
%    U=cat(3,-2,-7,-4);SL3(U),summary(U,1e-8)
%    missing parameter - default view:
%    1
%    1
%    1
%    CHM norm = 2.89807e-13
%    |1| norm = 4.57757e-16
%    defect   = 4
%    * this is a Butson BH(N=9, 18)!
%    #L       = 9
%
%
% only INTEGER triplets for odd values of d give 2U+CHM
% most examples work for |c|=2*k for ONLY particular values of k in N
% for example, not for c=12, cf. cat.m for other examples
%
% other findings:
% d=3 (N=9)  : (a,b,c) = (-5,  2, -2) ==> U in BH(9, 18),   defect(U)=4,   #L(U)=9
% d=5 (N=25) : (a,b,c) = (-6, -4, -2) ==> U in BH(25, 25),  defect(U)=16,  #L(U)=...(LONG CALCULATIONS)
% d=5 (N=25) : (a,b,c) = ( 4,  1,  2) ==> U in BH(25, 50),  defect(U)=16,  #L(U)=...(LONG CALCULATIONS)
% d=7 (N=49) : (a,b,c) = ( 1, -1,  2) ==> U in BH(49, 98),  defect(U)=36,  #L(U)=...(LONG CALCULATIONS)
% d=9 (N=49) : (a,b,c) = (-1, -2,  2) ==> U in BH(81, ...), defect(U)=136, #L(U)=...(LONG CALCULATIONS)
%
% for even values of d, cat map does not work:
% for d=4 (N=16), Z is not even close to zero..., best Z ~ 0.0618059
% for d=6 (N=26) : a=0.796323       b=-1.73002      c=-2,   Z=0.0300646     [RESET@#2215]
%
%                  a=1.00044        b=-2.50088      c=-2,   Z=0.0139092     [RESET@#26]
%                  >> cat(6, 1, -2.5, -2) is BH(36, 144) with defect=97 and "UNCOUNTABLE" #L
%
%                  a=-1.73972       b=-0.999434     c=-2,   Z=0.0121736     [RESET@#1]
%                  >> cat(6, -1, 3/4, -2) is BH(36, 288) with defect=97 and "UNCOUNTABLE" #L
%
%
% This ALWAYS converges to |c| approx. 2.
%
% ------------------------------------------------------------------------------

    addpath ../matrix_tool
    addpath ../CHM
    N = d^2;

    while 1
        a = a0;#randi(10)-10;
        b = b0;#randi(10)-10;
        c = c0;#randi(10)-10;
        mu = 1;
        Z_OPTIMAL = +Inf;
        m = mu0;

        while Z_OPTIMAL >1e-13 && m<1e+6
            m++;
            RESTORE_a = a;
            RESTORE_b = b;
            RESTORE_c = c;

            a += randn*mu; % a = -6;
            b += randn*mu;
            c += randn*mu;

            U = zeros(N, N);
            for j=1:N
            for k=1:N
                U(j, k) = (1/1)*exp((1j*pi/N)*(a*(j^2)+b*(k^2)+c*j*k));
            end, end

            Z = abs(SL(U) - 1) + abs(SL(R(U)) - 1) + abs(SL(T(U)) - 1) + nh(U);

            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z * 0.001;
                printf("a=%g \t b=%g \t c=%g,\t Z=%g \t [RESET@#%d]\n", a, b, c, Z, m);
                m = 0;
            else
                mu = rand;
                a = RESTORE_a;
                b = RESTORE_b;
                c = RESTORE_c;
            end
        end

        if Z_OPTIMAL < 1e-6
            fileName = strcat("cat_walk_", int2str(N), "_", datestr(now(), 30), ".data");
            save(fileName, "*");
        end

    end

end

