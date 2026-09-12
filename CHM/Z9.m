function Y = Z9(p)
% ------------------------------------------------------------------------------
% 2008-05-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
%
% ********************
% *                  *
% *  AUTHOR WANTED!  *
% *                  *
% ********************
%
% -----------------------------------------------------------------------------
% Generic values: d(Y) = 2 and #L(Y) = 58.
% For example, for p = 0.25 one recovers BH(9, 12) with d = 10 and #L = 12.
% -----------------------------------------------------------------------------

    if p(1) < 0.1 || p(1) > 0.33
        error("Wrong value of the free parameter!");
    end

    u2 = exp(2j*pi* p(1));

    p = 1-2*u2+u2^2-2*u2^3+u2^4+(-1+u2-u2^2)*sqrt((1+u2*(u2-3))*(1+u2+u2^2));
    u1 = p/2/u2;
    u4 = -1+u2*(p/2/u2 - p*(sqrt(p/u2)/sqrt(2)+1/sqrt(u2)+u2^(3/2))/(sqrt(p/u2)*u2/sqrt(2)-u2^(3/2)) + (sqrt(p/u2)/sqrt(2)+1/sqrt(u2)+u2^(3/2))^2);
    u3 = (sqrt(p/u2)/sqrt(2)+1/sqrt(u2)+u2^(3/2))/(sqrt(p/u2)/sqrt(2)-sqrt(u2));
    % it is not the only possible configuration for uj...
    % others are much more complicated...

    Y = [
        1                      1                      1                      1                      1            1                      1                      1                      1  ;
        1                     -1                 -u1*u2                     u2                     u1 -sqrt(u1*u2)                  u1*u2               u1*u2*u3                     u3' ;
        1              (u2*u3)^2 sqrt(u1)*u2^(3/2)*u3^2                 -u2*u3                u2*u3^2           u2                     -1             -u2^2*u3^2                     u3  ;
        1                     u2               u1*u2*u3                     -1                    u3'        u1*u2           -sqrt(u1*u2)                 -u1*u2                     u1  ;
        1             -(u2*u3)^2            -u1*u2*u3^2 sqrt(u1)*u2^(3/2)*u3^2 sqrt(u1)*u2^(3/2)*u3^2       -u1*u2               u1*u2*u3                     u4               u1*u2*u3  ;
        1                     u3 sqrt(u1)*u2^(3/2)*u3^2                u2*u3^2                     u2           u1                     u3'              u1*u2*u3             u1*u2*u3^2  ;
        1                u2*u3^2               u1*u2*u3                     u3             u1*u2*u3^2           u3'                    u1 sqrt(u1)*u2^(3/2)*u3^2                     u2  ;
        1                 -u2*u3             -(u2*u3)^2              (u2*u3)^2                     u3           -1                     u2 sqrt(u1)*u2^(3/2)*u3^2                u2*u3^2  ;
        1 sqrt(u1)*u2^(3/2)*u3^2                     u4             -(u2*u3)^2               u1*u2*u3     u1*u2*u3                 -u1*u2            -u1*u2*u3^2 sqrt(u1)*u2^(3/2)*u3^2  ;
    ];

end
