function Y = SR_9_search
% ------------------------------------------------------------------------------
% 2022-11-16 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2022-12-27
% 2023-01-02
%
% Searching for a self R-dual CHM of order 9 witch 54 independent parameters.
% For a general case with all 54 parameters the best solutions is around Z ~ 0.002...
% Simplified (symmetrized) version provides a numerical solution: X = R(X) in H(9).
%
% Slower convergence in comparison with "SG_9_search.m".
% ------------------------------------------------------------------------------

    np = 54;
    RESET_INDEX = 0;
    Z_OPTIMAL = +Inf;
    mu = 1;
    p = rand(1, np); % initial random phases
    k = 0;

    while (Z_OPTIMAL>5e-13) % && (k<1000000)
        RESTORE_p = p;
        p = mod(p + (randn(1, np) * mu), 1);
        X = exp(2j*pi*[
            p(1) p(2) p(3)      p(4) p(5) p(6)      p(7) p(8) p(9);
            p(4) p(5) p(6)      p(10) p(11) p(12)   p(13) p(14) p(15);
            p(7) p(8) p(9)      p(13) p(14) p(15)   p(16) p(17) p(18);

            p(19) p(20) p(21)   p(22) p(23) p(24)   p(25) p(26) p(27);
            p(22) p(23) p(24)   p(28) p(29) p(30)   p(31) p(32) p(33);
            p(25) p(26) p(27)   p(31) p(32) p(33)   p(34) p(35) p(36);

            p(37) p(38) p(39)   p(40) p(41) p(42)   p(43) p(44) p(45);
            p(40) p(41) p(42)   p(46) p(47) p(48)   p(49) p(50) p(51);
            p(43) p(44) p(45)   p(49) p(50) p(51)   p(52) p(53) p(54);
        ]);

        assert(norm(X - R(X), "fro") < 1e-12);

        Z = nh(X);
        if Z < Z_OPTIMAL
            printf("Z=%g \t mu=%g \t RESET@k=%d\n", Z, mu, k);
            k = RESET_INDEX;
            Z_OPTIMAL = Z;
            mu = Z * 0.001;
        else
            p = RESTORE_p;
        end
        k++;
    end % while ...

    Y = X;
    A = mod(angle(Y)+2*pi,2*pi)/2/pi;
    fileName = strcat("sRD_9_d", num2str(ud(Y,"S",1e-8)), "_", datestr(now(), 30), ".data");
    save(fileName, "A");

end

