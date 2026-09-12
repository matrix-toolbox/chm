function H = K6_3(p)
% ------------------------------------------------------------------------------
% 2009-08-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Family of CHM depending on 2 parameters found by B. Karlsson;
% [1] https://arxiv.org/abs/1003.4177
%
% Watch the scope of paramaters!
%
% >> H = K6_3(rand(1, 3));
% ------------------------------------------------------------------------------

     try
         theta = 2*pi*p(1);
         phi = p(2);
         psi = p(3);
     catch
         warning("No valid parameter provided! Continue with default [1 1 1] value!");
         theta = 1;
         phi = 1;
         psi = 1;
     end

     e = [1 0; 0  1];
     F2 = [1 1; 1 -1];
     Lambda = [
         cos(theta), exp(2i*pi* phi) * sin(theta);
         exp(-2i*pi * phi) * sin(theta), -cos(theta)
     ];
     A = F2 * (-e / 2 + 1i * sqrt(3) / 2 * Lambda);
     B = - F2 - A;
     z1 = exp(2i * pi* psi);

     if abs(abs(A(1, 2)) - abs(A(1, 1))) > 100 * eps && abs(abs(B(1, 2)) - abs(B(1, 1))) > 100 * eps
         Moebius         = ['(alpha * z^2 - beta) / (conj(beta) * z^2 - conj(alpha))'];
         Moebius_inverse = ['(conj(alpha) * z^2 - beta) ./ (conj(beta) * z^2 - alpha);'];
         alpha           = A(1, 2)^2;
         beta            = A(1, 1)^2;
         z               = z1;
         z3              = sqrt(eval(Moebius));
         alpha           = B(1, 2)^2;
         beta            = B(1, 1)^2;
         z4              = sqrt(eval(Moebius));
         z               = z3;
         z2              = sqrt(eval(Moebius_inverse));
         Z1              = [1 1; z1 -z1];
         Z2              = [1 1; z2 -z2];
         Z3              = [1 z3; 1 -z3];
         Z4              = [1 z4; 1 -z4];

         H = [
            F2 Z1        Z2       ;
            Z3 Z3*A*Z1/2 Z3*B*Z2/2;
            Z4 Z4*B*Z1/2 Z4*A*Z2/2;
         ];
     else
        error("Degenerated case not implemented!")
     end

end
