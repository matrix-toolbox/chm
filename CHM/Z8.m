function H = Z8(p)
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
% A 1-parameter family of CHM of order N = 8,
% with generic d(Y) = 5 and #L(Y) = 102.
%
% For p = 0.5 it intersects BH(8, 2) with d = 21 and #L = 2.
%
% >> H = Z8(rand);
% ------------------------------------------------------------------------------

  warning("Scope of the parameter is restricted!");


    if p(1) < 0.42063623618066176 || p(1) > 0.5793637638193383
        error("Wrong value of the free parameter!");
    end

    x1 = exp(2*pi*i * p(1));
    x4 = (-1-x1-x1^2-sqrt(1+2*x1+2*x1^2+2*x1^3+x1^4)) / x1;

    A1 = -7-8*x1-10*x1^2-10*x1^3-4*x1^4-x1^6;
    A2 = 17+33*x1+42*x1^2+47*x1^3+42*x1^4+21*x1^5+15*x1^6+6*x1^7+x1^9;
    A3 = 2-2*x1+13*x1^2+12*x1^3+25*x1^4+28*x1^5+25*x1^6+12*x1^7+13*x1^8-2*x1^9+2*x1^10;

    x7 = ( -2-2*x1-x1^3+ A1 / (A2+3*sqrt(3)*sqrt(-(1+x1)^2*( A3 ) ) )^(1/3) - ...
         (A2+3*sqrt(3)*sqrt(-(1+x1)^2*( A3) ) )^(1/3)) / 3;
    x6 = (-2*x7-2*x1*x7-x1^3*x7-x7^2-sqrt(-4*x1^3*x7+(2*x7+2*x1*x7+x1^3*x7+x7^2)^2)) / 2 / x7;
    x5 = - (x1*(1+x1^2+x1*x4)) / (x1^2+x1*x7+x6*x7);
    x3 = x7*x5 / x1;
    x2 = x6*x3/x1;

    H = [
        1  1         1          1         1         1    1         1  ;
        1 -1        x1'        x2'       x1        x3'  x4'       x5' ;
        1 x1         1   x1*x5/x2      x1^3  x1*x2/x3   x1  x1*x3/x5  ;
        1 x2  x2/x1/x5         -1  x1*x2/x3 -x1*x2/x3  -x2 -x2/x1/x5  ;
        1 x1'   1/x1^3   x3/x1/x2         1  x5/x1/x3   x1' x2/x1/x5  ;
        1 x3  x3/x1/x2  -x3/x1/x2  x1*x3/x5        -1  -x3 -x1*x3/x5  ;
        1 x4        x1'       -x2'       x1       -x3'   1       -x5' ;
        1 x5  x5/x1/x3  -x1*x5/x2  x1*x5/x2 -x5/x1/x3  -x5        -1  ;
    ];

end
