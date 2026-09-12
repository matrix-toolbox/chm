function Y = T6(p)
% ------------------------------------------------------------------------------
% 2008-10-20 WB
% 2021-05-15 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Non-affine family of CHM of order N = 6.
%
% Parameter m (or p) is a root of palindromic minimal polynomial:
% x^12 + 3x^11 + 8x^10 + 13x^8 - 3x^7 + 20x^6 - 3x^5 + 13x^4 + 8x^2 + 3x + 1
%  * sum of roots = -3
%  * product of roots = 1
%
% Generic values: d(Y) = 4, #L(Y) = 22
% p = 1/4 --> d = 4, #L = 4, symmetric BH(6, 4)
% p = 1/2 --> d = 4, #L = 4, Hermitian BH(6, 4)
% p = 3/4 --> d = 4, #L = 4, symmetric BH(6, 4)
%
% >> Y = T6(rand/2 + 1/4);
% ------------------------------------------------------------------------------

    if p(1) < 1/4 || p(1) > 3/4
        error("Scope of the parameter does not cover entire circle; p in [1/4, 3/4].");
    end

    a = exp(2j * pi * p);
    % m = (- 2 - a' - a - sqrt(a' * a' + a * a + 2 + 4 * (a' + a))) / (2 * a');
    % p = (- 2 - a' - a + sqrt(a' * a' + a * a + 2 + 4 * (a' + a))) / (2 * a');
    b = a^2; % b = conj(- 1 - 2 * a' - c' - d');
    c = (-2*a - 1 - b - sqrt(1+a^4+2*a^2+4*a+4*a^3))/2;
    d = (-2*a - 1 - b + sqrt(1+a^4+2*a^2+4*a+4*a^3))/2;

    % q = 7 (11 ...)
    % (c*d)^q = 1
    % c*d = exp(12*1j*pi/q)

    Y = [
        1 1 1  1  1  1;
        1 a b  c  a  d;
        1 b a  a  c  d;
        1 d a -a -1 -d;
        1 a d -1 -a -d;
        1 c c -c -c -1;
    ];

end

% ------------------------------------------------------------------------------
% MATHEMATICA version:
%
% (* scope of parameter $p$ in [1/4, 3/4] *)
% Clear[a, c, d, p];
% a[p_] := Exp[2 Pi I p];
% c[p_] := (-2 a[p] - 1 - a[p]^2 - 
%      Sqrt[1 + a[p]^4 + 2 a[p]^2 + 4 a[p] + 4 a[p]^3])/2;
% d[p_] := (-2 a[p] - 1 - a[p]^2 + 
%      Sqrt[1 + a[p]^4 + 2 a[p]^2 + 4 a[p] + 4 a[p]^3])/2;
% Plot[If[Abs[Abs@c[p] - 1] < 0.0001, 1, 0.0], {p, 0, 1}, 
%  PlotStyle -> {Directive[Blue, Thick]}, Filling -> Axis, 
%  FillingStyle -> LightBlue, PerformanceGoal -> "Quality", 
%  AspectRatio -> 1/GoldenRatio,
%  TicksStyle -> Directive[Black, Bold, 18], 
%  AxesLabel -> {Style["p", Bold, 18], Style[Abs[c], Bold, 18]}, 
%  Ticks -> {{0, 0.25, 0.5, 0.75, 1}, {0, 1}}]
%
% T6[p_] := {
%   {1, 1, 1, 1, 1, 1},
%   {1, a[p], a[p]^2, c[p], a[p], d[p]},
%   {1, a[p]^2, a[p], a[p], c[p], d[p]},
%   {1, d[p], a[p], -a[p], -1, -d[p]},
%   {1, a[p], d[p], -1, -a[p], -d[p]},
%   {1, c[p], c[p], -c[p], -c[p], -1}
%   };
%
% Y = T6[RandomReal[{0.25, 0.75}]];
% Y.ConjugateTranspose[Y] // MatrixForm // N // Chop
% UnitaryMatrixQ[Y/Sqrt[6]]
% Abs@Y // MatrixForm // N
% UnitaryDefect[Y]
% ------------------------------------------------------------------------------
