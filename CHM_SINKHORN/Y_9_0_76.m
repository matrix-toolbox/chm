function Y = Y_9_0_76
% ------------------------------------------------------------------------------
% 2023-01-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% An Isolated and symmetric CHM matrix of order 9 with #L(Y) = 76
% found by the Sinkhorn algorithm.
%
% [1] https://arxiv.org/abs/2204.11727
% ------------------------------------------------------------------------------

    printf("Wait, calculating solution...");
    do
        [u info]=fsolve(@uc, exp(2j*pi*rand(1, 3)));
        for j=1:64
            [u info]=fsolve(@uc, u);
        end
        s = abs(abs(u(1)) - 1) + abs(abs(u(2)) - 1) + abs(abs(u(3)) - 1);
        printf(".");
    until s < 1e-13
    printf(" Solved!\n");

    x = u(1)
    y = u(2)
    z = u(3)

    a = 2*(y + z) + (2 + x/y)*x/z + (2*x^2)/((y + z)*(x + y*z)) - 1;
    b = -(x/y) - y - x/z - z;
    c = 2*x + ((2*x)/(y*z) + 3)*x^2/(y*z) - 1;

    Y = [
        1, 1, 1, 1, 1, 1, 1, 1, 1;
        1, a, b, b, x/y, x/y, -(x^2/(y*z)), x^2/(b*y*z), x^2/(b*y*z);
        1, b, x, y, x/y, z, -x, -1, x/z;
        1, b, y, x, z, x/y, -x, x/z, -1;
        1, x/y, x/y, z, x^2/y^2, -(z/y), -(x^2/(y*z)), x/y, x^2/(y^2*z);
        1, x/y, z, x/y, -(z/y), x^2/y^2, -(x^2/(y*z)), x^2/(y^2*z), x/y;
        1, -(x^2/(y*z)), -x, -x, -(x^2/(y*z)), -(x^2/(y*z)), c, -(x^3/(y^2*z^2)), -(x^3/(y^2*z^2));
        1, x^2/(b*y*z), -1, x/z, x/y, x^2/(y^2*z), -(x^3/(y^2*z^2)), x^3/(y^2*z^2), x^2/(y*z^2);
        1, x^2/(b*y*z), x/z, -1, x^2/(y^2*z), x/y, -(x^3/(y^2*z^2)), x^2/(y*z^2), x^3/(y^2*z^2);
    ];


end



function u=uc(v)
	u(1) = v(1) - (3*v(2)*v(3))/v(1) + ((1 - 2*v(1))*v(2)^2*v(3)^2)/v(1)^3 - 2;
	u(2) = 1/v(2) + v(1)/(v(2)*v(3)^2) + 1/v(3) + v(1)/(v(2)^2*v(3)) - v(1)/((v(2) + v(3))*(v(1) + v(2)*v(3)));
	u(3) = 3 + (v(1)*(1 - v(2) + v(3)))/(v(2)*v(3)) + (v(2) + (-1 + v(2))*v(3))/v(1);
end
