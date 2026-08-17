function Y = Y_10_4_114
% -----------------------------------------------------------------------------
% 2022-08-23 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% A 4-parametric non-affine family of CHM of order N = 10 found by the Sinkhorn algorithm.
% It has generic defect = 4 and #L = 114.
%
% [1] https://arxiv.org/abs/2204.11727
% -----------------------------------------------------------------------------

    #    a =  0.980248823743152 - 0.197768156056959i;
    #    b =  0.871151522932081 + 0.491014280946202i;
    #    c =  0.993007570197494 + 0.118050690512443i;
    #    d =  9.094145038285016e-02 - 9.958562409315227e-01i;

    # original values: (?)
    # a = -0.9988874350068101 - 0.04715816138820191 I;
    # b = -0.9990958699863974 + 0.04251402799221279 I;
    # c = 0.5775433708165421 + 0.8163600032007126 I; # p3 = 0.95508;
    # d = -0.4224350177633097 - 0.9063932125558496 I;

    #a=exp(2j*pi/20); % 20
    #b=exp(2j*pi/18);

    w=exp(2j*pi/6);


    Z_OPTIMAL = +Inf;
    Z = +Inf;
    mu = 1;

    p1 = rand;
    p2 = rand;
    p3 = rand;
    p4 = rand;


    while Z>5e-14
        RESTORE_p1 = p1;
        RESTORE_p2 = p2;
        RESTORE_p3 = p3;
        RESTORE_p4 = p4;

        p1 = mod(p1 + randn*mu, 1);
        p2 = mod(p2 + randn*mu, 1);
        p3 = mod(p3 + randn*mu, 1);
        p4 = mod(p4 + randn*mu, 1);

        a = exp(2j*pi*p1);
        b = exp(2j*pi*p2);
        c = exp(2j*pi*p3);
        c = exp(1j*p3);
        d = exp(2j*pi*p4);

        Y=[
           1           1            1           1           1            1            1             1             1           1      ;
           1           1            a          -a           w/a         -w/a          w^2           w^2          -w          -w      ;
           1           b           -a*b         a          -1           -b            a*b*c*w^2    -a*b*c*w^2     a*b        -a      ;
           1          -b            b          -1          -w/a         -w*b/a        w*b*c        -w*b*c         w*b/a       w/a    ;
           1           w/b         -1          -w/b        -w^2/a/b      w/a          c            -c             w^2/a/b    -w/a    ;
           1          -w/b         -a          -w*a/b       w/b         -1           -w*a*c         w*a*c         w*a/b       a      ;
           1           w^2         -a*b*d*w^2  -w*a*d      -d            w*b*d        w^5           w^5          -w           w^2    ;
           1           w^2          a*b*d*w^2   w*a*d       d           -w*b*d        w^5           w^2          -w           w^2    ;
           1          -w            a*b         w*a/b       w^2/a/b      w*b/a       -w            -w             w^2        -w      ;
           1          -w           -b           w/b        -w/b          b            w^2           w^2          -w           1      ;
        ];

        Z = nh(Y);
        if Z < Z_OPTIMAL
            Z_OPTIMAL = Z;
            mu = Z_OPTIMAL * 0.01;
        else
            p1 = RESTORE_p1;
            p2 = RESTORE_p2;
            p3 = RESTORE_p3;
            p4 = RESTORE_p4;
        end


    end


    p1
    p2
    p3
    p4

    a,b,c,d

end

