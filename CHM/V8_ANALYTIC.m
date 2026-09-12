function V8_ANALYTIC
% ------------------------------------------------------------------------------
% 2021-09-11 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% One possible configuration for analytic values of parameters a, b and c
% from V_8^{(0)}.
%
% Based on notes of Kulkarni (note the mistake in one formula!)
% https://scholarworks.umt.edu/cgi/viewcontent.cgi?article=1071&context=tme
% ------------------------------------------------------------------------------

    format long g
    [z convergence_info] = fsolve(@get_solution, exp(2j*pi*rand));
    % z
    % convergence_info

    % choose appropriate vector for polynomial coefficients in function "get_solution" below

    v1 = -z/2 -sqrt(z^2-4)/2
    v2 = +z/2 -sqrt(z^2-4)/2

    % possible solutions (Mathematica can provide complicated but analytical form)
    % for the 1st vector:
    % v1 =     0.7799568409948809 -     0.6258333054233315i
    % v2 =    -0.7799568410012712 -     0.6258333054284589i
    %
    % v1 =    -0.1941704031915285 -     0.9809678155375114i
    % v2 =     0.1941704030758034 -     0.9809678149528567i
    %
    % for the 2nd vector:
    % v1 =   -0.03255228911764775 +     0.9994700339262065i
    % v2 =    0.03255228910972936 +     0.9994700336830841i
    %
    % v1 =     0.6183387234035549 -     0.7859117055500533i
    % v2 =    -0.6183387321957556 -     0.7859117167249856i
    %
    % then, for {a, b, c} = {-v1/v4, -v1, v3}
    % one obtains analytic form of
    % V = {
    % {-1, -1,  b,  b,  c,  c,  a,  a},
    % {-1,  b, -1,  c,  b,  a,  c, -a},
    % { b, -1,  c, -1,  a,  b, -a,  c},
    % { b,  c, -1,  a, -1, -a,  b, -c},
    % { c,  b,  a, -1, -a, -1, -c,  b},
    % { c,  a,  b, -a, -1, -c, -1, -b},
    % { a,  c, -a,  b, -c, -1, -b, -1},
    % { a, -a,  c, -c,  b, -b, -1,  1}
    % };

    % 2021-09-13
    % formulas for eta, alpha_j and beta_j can be simplified further...
    % eta = 464 - sqrt(512)/2;
    % alpha1 = -4+sqrt(eta);
    % alpha2 = -16+4*sqrt(eta)-1312/sqrt(eta);
    % alpha3 = -228+eta/2+2*sqrt(eta)*(eta-440)/(328-eta);
    % beta1 = -4-sqrt(eta)
    % beta2 = -16-4*sqrt(eta)+1312/sqrt(eta)
    % beta3 = -228+eta/2-2*sqrt(eta)*(eta-440)/(328-eta)

end

% ------------------------------------------------------------------------------

function y = get_solution(x)
    a0 = -16;
    a1 =  256;
    a2 = -96;
    a3 = -896;
    a4 = -696;
    a5 = -96;
    a6 =  56;
    a7 =  16;

    F2  = (4*a6-a7^2)/8;
    F1  = (a5-a7*F2)/2;
    F0  = (a4-a7*F1-F2^2)/2;
    F3  = (a3-a7*F0-2*F1*F2)/2;
    F4  = a2 - F1^2 - 2*F0*F2;
    F5  = (a1-2*F0*F1)/2;
    F6  = a0-F0^2;
    F7  = (16*F2-a7^2)/16;
    F8  = (a7*F3-2*F4)/2;
    F9  = (a7*F8+8*F5-4*F3*F7)/(a7*F7-4*F1);
    F10 = -(a7*F3^2+4*F3*F8)/(a7*F7-4*F1);
    F11 =  4*F3^3/(a7*F7-4*F1);
    F12 =  (4*a7^2*F0-16*F1^2-8*a7*F3)/(a7^2);
    F13 =  (16*F3^2+64*F1*F5-4*a7^2*F6-32*a7*F0*F3)/a7^2;
    F14 =  (32*a7*F3*F6+64*F0*F3^2-64*F5^2)/a7^2;
    F15 = -64*F3^2*F6/a7^2;

    hh = (F10+F9*F12-F13-F9^2);
    h0 = (F11*F12-F15-F9*F11)/hh;
    h1 = (F11+F10*F12-F14-F9*F10)/hh;

    g = -(h1+sqrt(h1^2-4*h0))/2;

    alpha1 = F2+sqrt(g)
    alpha2 = F1+a7*sqrt(g)/4-F3/sqrt(g)
    alpha3 = F0+g/2+sqrt(g)*(8*F5-4*F1*g)/(8*F3-2*a7*g)

    beta1 = F2-sqrt(g)
    beta2 = F1-a7*sqrt(g)/4+F3/sqrt(g)
    beta3 = F0+g/2-sqrt(g)*(8*F5-4*F1*g)/(8*F3-2*a7*g)

    % F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15
    % a7/2, alpha1, alpha2, alpha3

    % choose either vector to get different pairs of solutions
%    c = [a7/2, alpha1, alpha2, alpha3]
    c = [a7/2, beta1, beta2, beta3]
    y = x^4 + c(1)*x^3 + c(2)*x^2 + c(3)*x + c(4);
end

% ------------------------------------------------------------------------------
% 2020-07-05
% ALTERNATIVE ANALYTIC APPROACH
%
% https://dl.acm.org/doi/10.1145/365758.365809
% http://pari.math.u-bordeaux.fr/dochtml/html-stable/Arithmetic_functions.html#factormod
% https://zims-en.kiwix.campusafrica.gos.orange.com/wikipedia_en_all_nopic/A/Octic_equation
% https://ask.sagemath.org/question/39772/factoring-polynomials-with-symbolic-expressions/
% https://cocalc.com/projects/44689845-cc82-434c-8c75-d50017c39456/files/Welcome%20to%20CoCalc.ipynb?session=default
% 
% @interact
% def _(p = input_box(default=x^8 + 16*x^7 + 56*x^6 - 96*x^5 - 696*x^4 - 896*x^3-96*x^2+256*x-16, label = 'P(x) = ')):
%     q.<x> = PolynomialRing(QQ, 'x')
%     r.<x> = PolynomialRing(RR, 'x')
%     c.<x> = PolynomialRing(CC, 'x')
%     a.<x> = PolynomialRing(AA, 'x')
%     qp = PolynomialRing(QQ, 'x')(p)
% #    print(factor(qp))
% #    rp = r(p)
% #    print(factor(rp))
% #    cp = c(p)
% #    print(factor(cp))
%     ap = a(p)
%      print(factor(ap))
% 
%     f=factor(ap)
%     c = list(f)[0][0].coefficients()[0]
%     print(c)
%     bb = c.radical_expression()
%     print(bb)
%      print(bb.degree())
% ------------------------------------------------------------------------------
