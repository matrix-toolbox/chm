function Y = Y_8_5_82
% -----------------------------------------------------------------------------
% 2022-08-27 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% We cannot say that parameters {a,b,c,d,e} are fully independent, however...
% Generic defect(Y) = 5 and #L(Y) = 82.
% It is probably a sub-orbit of F8.
% -----------------------------------------------------------------------------

    aa = 2*pi*rand;
    ab = 2*pi*rand;
    ac = 2*pi*rand;
    ad = 2*pi*rand;
    ae = 2*pi*rand;
    w = exp(2j*pi/12);

    Z_OPTIMAL = +Inf;
    Z = +Inf;
    mu = 1;

    while Z>5e-14
        RESTORE_aa = aa;
        RESTORE_ab = ab;
        RESTORE_ac = ac;
        RESTORE_ad = ad;
        RESTORE_ae = ae;

        aa = mod(aa + randn*mu, 1);
        ab = mod(ab + randn*mu, 1);
        ac = mod(ac + randn*mu, 1);
        ad = mod(ad + randn*mu, 1);
        ae = mod(ae + randn*mu, 1);
        a = exp(2j*pi*aa);
        b = exp(2j*pi*ab);
        c = exp(2j*pi*ac);
        d = exp(2j*pi*ad);
        e = exp(2j*pi*ae);

        Y = [
            1       1       1       1       1       1       1       1        ;
            1       1       1       1      -1      -1      -1      -1        ;
            1       1      -1      -1       b       b      -b      -b        ;
            1       1      -1      -1      -b      -b       b       b        ;
            1      -1       c      -c       d      -d       e      -e        ;
            1      -1       c      -c      -d       d      -e       e        ;
            1      -1      -c       c       a      -a       w      -w        ;
            1      -1      -c       c      -a       a      -w       w        ;
        ];

        Z = nh(Y);
        if Z < Z_OPTIMAL
            Z_OPTIMAL = Z;
            mu = Z_OPTIMAL * 0.01;
        else
            aa = RESTORE_aa;
            ab = RESTORE_ab;
            ac = RESTORE_ac;
            ad = RESTORE_ad;
            ae = RESTORE_ae;
        end


    end


    a
    b
    c
    d
    e

end



