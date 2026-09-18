function [Y, x, y, z] = BH36_6_D_2U_NA(p, t)
% ------------------------------------------------------------------------------
% 2026-08-23
% Wojciech Bruzda | name@cft.edu.pl : name = w.bruzda
% ------------------------------------------------------------------------------
% one-parameter non-affine extension of the Dita family
%
% >> [Y, x, y, z] = BH36_6_D_2U_NA(p, t)
%
% returns Y = BH36_6_D_2U_A(p) .* Z(x, y, z) : x = exp (1i*t)
%
% which is unimodular and two-unitary for every real p (15 entries, see BH36_6_D_2U_A) and every t in the open interval
%
% (t_minus, t_plus) = (-4.4865188906..., 0.2977286858...)
%
% t_plus  = pi / 3 - acos(sqrt(3) - 1)
% t_minus = pi / 3 + acos(sqrt(3) - 1) - 2 * pi
%
% t = 0 gives x = y = z = 1, the multiplier is all ones and the affine family is recovered
% for generic t the phases are not roots of unity, so Y is a two-unitary complex Hadamard matrix that is NOT of Butson type:
% this is the direction that no Dita deformation can produce
%
% Z depends on the entries only through:
% m = mod(a + c, 6)
% nn = mod(b + d, 2)
% Z(6 * a + b + 1, 6 * c + d + 1) = Phi(m + 1, nn + 1)
% Phi = [y 1; 1 z; 1 x; x z; x y; y z]
%
% note that nn sees b and d only through their parities, i.e., only through the Z_2 part of the Chinese remainder splitting used in BH36_6_D_2U
% the deformation acts on the twelve-dimensional core and leaves the order-three Fourier blocks untouched
% it couples the legs A and C, which the Dita arrays psi(J, be) never do -- hence non-affine!
%
% z = z(x) is the analytic root through z(1) = 1 of F(x, z) = A(x) * z^2 + B(x) * z + C(x) = 0
% where A(x) = 1 + w + (w - 2) * x
%       B(x) = x^2 - 4 * w * x + w^2
%       C(x) = x * ((2 * w-1) * x + 2 - w) : w = exp(1j * pi / 3)
%
% and y = w * x * z / (x * (z + 2 * w - 1) - w * z)
%
% WARNING
% do not evaluate the closed form with the principal square root
% arg(Delta(exp(1j * t))) crosses the branch cut inside the admissible interval, which silently swaps the two roots of F
% the discrepancy against the correct branch reaches 1.15
% hence, the square root is continued analytically along the arc from 0 to t below, normalised by sqrt (Delta(1)) = 2 - w
% and the result is checked against F(x, z) = 0
%
% >> Y = BH36_6_2U_D_NA(zeros(1, 15), 0.2);
% ------------------------------------------------------------------------------

    if nargin < 1 || isempty (p), p = zeros (15, 1); end
    if nargin < 2 || isempty (t), t = 0; end
    if ~isreal (t) || ~isscalar (t)
        error ("paramater t must be a real scalar!");
    end

    [x, y, z] = na_branch(t);

    Phi = [y 1; 1 z; 1 x; x z; x y; y z];
    Z = zeros (36, 36);
    for a = 0:5
    for b = 0:5
    for c = 0:5
    for d = 0:5
        Z(6*a+b+1, 6*c+d+1) = Phi(mod (a+c, 6) + 1, mod (b+d, 2) + 1);
    end, end, end, end

     Y = BH36_6_D_2U_A(p) .* Z;
end

function [x, y, z] = na_branch(t)
% analytic branch (x,y,z) = (exp(1i * t), y(x), z(x)) through (1, 1, 1)
    w  = exp (1j * pi / 3);
    A  = @(x) 1 + w + (w - 2) * x;
    B  = @(x) x.^2 - 4 * w * x + w^2;
    C  = @(x) x.*((2 * w - 1) * x + 2 - w);
    DD = @(x) B(x).^2 - 4 * A(x) .* C(x);

    tp = pi / 3 - acos (sqrt (3) - 1);
    tm = pi / 3 + acos (sqrt (3) - 1) - 2 * pi;
    if t > tp || t < tm
        warning("outside", "t = %.6f lies outside (%.6f, %.6f): |z| ~= 1 and Y is not Hadamard", t, tm, tp);
    end

    x = exp (1j * t);
    if t == 0
        y = 1; z = 1; return;
    end

    N  = max(512, ceil(2000 * abs(t))); % analytic continuation of sqrt(Delta)
    th = linspace(0, t, N);
    Dv = DD (exp(1j * th));
    ph = unwrap(angle(Dv));
    s  = sqrt(abs(Dv(end))) * exp(0.5j * ph(end));
    z = (-B(x) - s) / (2 * A(x));
    d = x * (z + 2 * w - 1) - w * z;
    y = w * x * z / d;

    if abs (A(x) * z^2 + B(x) * z + C(x)) > 1e-10
        error ("branch continuation failed, |F| = %.3e", abs (A(x) * z^2 + B(x) * z + C(x)));
    end
end
