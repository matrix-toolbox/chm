function Y = T8(p)
% ------------------------------------------------------------------------------
% 2018-01-06 WB
% 2018-02-24 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 1-parametric non-affine family of CHM of order N = 8;
% * generic defect=3
% * p = 1/2, 1, 13/10  =>  defect=7  and T in BH(8, 20)
% * p = 4/5            =>  defect=11 and T in BH(8, 20)
%
% >> Y = T8_1(rand)
% ------------------------------------------------------------------------------

  g1 = 0.38629047387009;           % check: p = phase of x
  g2 = 0.8301134750001507;
  g3 = 0.9700212764097545;
  g4 = 1.413709525400110;
  if ~((p > g1 && p < g2) || (p > g3 && p < g4))
     warning("Wrong parameter! Continue with random (and valid) value...")
     u1 = rand;
     u2 = rand;
     u3 = fix(rand*2);
     r1 = u1 * g2 - (u1 - 1) * g1; % 1st interval
     r2 = u2 * g4 - (u2 - 1) * g3; % 2nd interval
     p = u3 * r1 + (1 - u3) * r2;
  end

  x = exp(i*pi*p);                 % free parameter

  s5=sqrt(5);                      % constants
  s2=sqrt(2);
  zeta1 = -i*sqrt(10-2*s5);
  zeta2 = +i*sqrt(10-2*s5);
  zeta4 = +i*sqrt(10+2*s5);
  
  x1 = 5 - 3*s5 + 4*sqrt(5-2*s5) + (1+i)/x + zeta2;
  x2 = 45-9 * i-(17-7*i)*s5 - 2*sqrt(620-500*i-(285-234*i)*s5) - (1+i)*x1/x;
  x3 = 10*i - 10 + (6+-2*i)*s5 + 4*sqrt(10-10*i-(5-4*i)*s5) + (1+i)*x2/x;
  x4 = s5 - 1 + i*s2*sqrt(5+s5) + (1+i)*x3/x;
  x5 = (1+i)*sqrt(1-s5+zeta4) - (1-i)*(2*x^2+sqrt(x4)') - x*(1-4*i+s5+zeta1);
  x6 = s5 - 1 - i - sqrt(5-2*s5) + (1-i)*x;
  y = x5 / x6 / 4;
  z = x * ((1+i)*exp(i*pi*4/5) - x - y - 1 + i) / (x-y);


%  q = (1/2)*(exp(i*pi/5) + exp(i*pi*7/10) - 1 + i)*(1/x-1/y)^(-1) + ...
%      (1/2)*(1/x-1/y)^(-2)*sqrt((1/x-1/y)^2*(i-1)*(2*exp(i*pi/5)+2*exp(i*pi*7/10) + ...
%      (1+s5+4*i+zeta2)*(1/x+1/y) + (1+i)*(2*(1/x+1/y)^2-exp(i*pi*9/10)+i)));
  q = -i*x/z;
  u = (q/x) * (z*x-x*y-y*z-x^2)/(1+i-exp(i*pi*3/10)-exp(i*pi*4/5));

  assert(abs(abs(y) - 1.0) < 1e-13)
  assert(abs(abs(z) - 1.0) < 1e-13)
  assert(abs(abs(u) - 1.0) < 1e-13)

  T = [
    0 0 0 0 0  0  0  0;
    0 0 0 0 0  3 15 18;
    0 8 8 8 8 13 15  8;
    0 2 2 2 2 17  5  2;
    0 0 0 0 0  7  5 12;
    0 0 0 0 0  0 10 10;
    0 0 0 0 0 10  0 10;
    0 0 0 0 0 10 10  0;
  ] / 10;

  Y = exp(i * pi * T) .* [
    1    1    1           1           1 1 1 1;
    1    x    y           z      -y*z/x 1 1 1;
    1 -u/y  u/x u*x/q/q/y/z     u/q/q/z 1 1 1;
    1    y    x      -y*z/x           z 1 1 1;
    1  u/x -u/y     u/q/q/z u*x/q/q/y/z 1 1 1;
    1    q   -q        -1/q         1/q 1 1 1;
    1    u   -u       u/q/q      -u/q/q 1 1 1;
    1  u/q  u/q        -u/q        -u/q 1 1 1;
  ];

end


% ------------------------------------------------------------------------------
% MATHEMATICA version:
% // sorry Stephen, you did a marvelous job inventing the language of Mathematica, but you totally f*cked up its (GUI) interface...
% 
% 
% Sqrt[1 - (7 + 2 Sqrt[5])/29.0] + I Sqrt[(7 + 2 Sqrt[5])/29.0];
% ArcTan[Sqrt[(7 + 2 Sqrt[5])/29]/Sqrt[1 - (7 + 2 Sqrt[5])/29]]/2/Pi // FullSimplify;
% a = -(ArcCot[1 - Sqrt[5]]/(2 \[Pi]));
% X = 1/4 Sqrt[8 - Sqrt[50 - 6 Sqrt[5]]] + I 1/4 Sqrt[8 + Sqrt[50 - 6 Sqrt[5]]];
% Z = ((-1 + I) - (-1)^(3/10) + (-1)^(4/5) - X + I Conjugate[X])/(1 - I X^2);
% W = -I Z X^2;
% T = ( {
%     {1, 1, 1, 1, 1, 1, 1, 1},
%     {1, W, Z, X, -I/X, Exp[Pi I 3/10], -I, Exp[Pi I 9/5]},
%     {1, -X, I/X, -I/Z, 1/W Exp[(3 Pi I)/2], Exp[Pi I 13/10], -I, Exp[Pi I 4/5]},
%     {1, Z Exp[(I Pi)/5], W Exp[(I Pi)/5], -(I/X) Exp[(I Pi)/5], X Exp[(I Pi)/5], Exp[Pi I 17/10], +I, Exp[Pi I 1/5]},
%     {1, -(I/X) Exp[(I Pi)/5], X Exp[(I Pi)/5], I/W Exp[(I Pi)/5], -(I/Z) Exp[(6 Pi I)/5], Exp[Pi I 7/10], +I, Exp[Pi I 6/5]},
%     {1, Exp[2 Pi I a], Exp[2 Pi I (a + 1/2)], Exp[2 Pi I (1/2 - a)], Exp[2 Pi I (1 - a)], +1, -1, -1},
%     {1, Exp[2 Pi I (a + 1/10)], Exp[2 Pi I (a + 3/5)], Exp[2 Pi I (11/10 - a)], Exp[2 Pi I (3/5 - a)], -1, +1, -1},
%     {1, Exp[Pi I 1/5], Exp[Pi I 1/5], Exp[Pi I 6/5], Exp[Pi I 6/5], -1, -1, +1}
%    } );
% 
% 
% \[Alpha] = Arg[Sqrt[(22 - 2 Sqrt[5])/29] + I Sqrt[(7 + 2 Sqrt[5])/29]] // FullSimplify; (* - ArcCot[1 - Sqrt[5]]*)
% Subscript[c, 1] = 1/4 Sqrt[8 - Sqrt[50 - 6 Sqrt[5]]] + I 1/4 Sqrt[8 + Sqrt[50 - 6 Sqrt[5]]];
% Subscript[c, 2] = ((-1 + I) - (-1)^(3/10) + (-1)^(4/5) - Subscript[c, 1] + I/Subscript[c, 1])/(1 - I Subscript[c, 1]^2);
% Subscript[c, 3] = -I Subscript[c, 2] Subscript[c, 1]^2;
% T = ( {
%     {0, 0, 0, 0, 0, 0, 0, 0},
%     {0, 0, 0, 0, 0, 3, 15, 18},
%     {0, 0, 0, 0, 0, 13, 15, 8},
%     {0, 2, 2, 2, 2, 17, 5, 2},
%     {0, 2, 2, 2, 2, 7, 5, 12},
%     {0, 0, 0, 0, 0, 0, 10, 10},
%     {0, 2, 2, 12, 12, 10, 0, 10},
%     {0, 2, 2, 12, 12, 10, 10, 0}
%    } );
% Subscript[T, c] = ( {
%     {1, 1, 1, 1, 1, 1, 1, 1},
%     {1, Subscript[c, 3], Subscript[c, 2], Subscript[c, 1], -(I/Subscript[c, 1]), 1, 1, 1},
%     {1, -Subscript[c, 1], I/Subscript[c, 1], -(I/Subscript[c, 2]), -(I/Subscript[c, 3]) , 1, 1, 1},
%     {1, Subscript[c, 2], Subscript[c, 3], -(I/Subscript[c, 1]), Subscript[c, 1] , 1, 1, 1},
%     {1, -(I/Subscript[c, 1]), Subscript[c, 1], I/Subscript[c, 3], I/Subscript[c, 2], 1, 1, 1},
%     {1, E^(I \[Alpha]), -E^(I \[Alpha]), -E^(-I \[Alpha]), E^(-I \[Alpha]), 1, 1, 1},
%     {1, E^(I \[Alpha]), -E^(I \[Alpha]), -E^(-I \[Alpha]), E^(-I \[Alpha]), 1, 1, 1},
%     {1, 1, 1, 1, 1, 1, 1, 1}
%    } );
% Subscript[T, 8] = Subscript[T, c]*Exp[I Pi T 1/10];
% Subscript[T, 8] // MatrixForm // N
% Subscript[T, 8].ConjugateTranspose[Subscript[T, 8]] // MatrixForm // N // Chop
% 
% 
% (* another candidate expressed in the analytical form *)
% (* with fixed unities at Subscript[T, 7,2] and Subscript[T, 7,3] *)
% x = 1/2 (1 - Sqrt[5] + Sqrt[2 (-2 + Sqrt[5])]) - I Sqrt[1/2 + Sqrt[-11 + 5 Sqrt[5]]];
% (*t=(-1+(-1)^(7/10)+(1+\[ImaginaryI]) (1+(-1)^(7/10)) x)/((-\[ImaginaryI]+(-1)^(1/5)) x);*)
% t = I - 1 - Cot[3/20 Pi] 1/x;
% (*z=((1+\[ImaginaryI]) (\[ImaginaryI]+(-1)^(4/5))+(-1+(-1)^(3/10)) x)/(1+(-1)^(3/10));*)
% z = I - 1 + I Tan[3/20 Pi] x;
% (*a=1-Sqrt[5]+Sqrt[5-2 Sqrt[5]]-(1-\[ImaginaryI])x;*)
% a = -Tan[3/20 Pi] - (1 - I) x;
% T690 = ( {
%     {1, 1, 1, 1, 1, 1, 1, 1},
%     {1, +x, I x Exp[I Pi 4/5], z, -I z Exp[I Pi 4/5], Exp[I Pi 3/10], -I, Exp[I Pi 9/5]},
%     {1, -(I/x), -(1/x) Exp[I Pi 4/5], t, +I t Exp[I Pi 4/5], Exp[I Pi 13/10], -I, Exp[I Pi 4/5]},
%     {1, -I x, x Exp[I Pi 1/5], z Exp[I Pi 1/2], -I z Exp[I Pi 7/10], Exp[I Pi 17/10], +I, Exp[I Pi 1/5]},
%     {1, -(1/x), I/x Exp[I Pi 1/5], I t, -t Exp[I Pi 1/5], Exp[I Pi 7/10], +I, Exp[I Pi 6/5]},
%     {1, +a, -a, -(1/a), +(1/a), +1, -1, -1},
%     {1, -1, +1, -(1/a^2), +(1/a^2), -1, +1, -1},
%     {1, -(1/a), -(1/a), +(1/a), +(1/a), -1, -1, +1}
%    } );
% 
% 
% (* matrix Subscript[T, 8]^(1)(x) *)
% (* 2018-02-24 W. Bruzda (improved by anonymous referee) *)
% 
% 
% (* constants *)
% Subscript[\[Zeta], 1] = Root[80 + 20 #1^2 + #1^4 &, 1];
% Subscript[\[Zeta], 2] = Root[80 + 20 #1^2 + #1^4 &, 2];
% Subscript[\[Zeta], 3] = Root[80 + 20 #1^2 + #1^4 &, 3];
% Subscript[\[Zeta], 4] = Root[80 + 20 #1^2 + #1^4 &, 4];
% 
% 
% (* draw free parameter at random *)
% (* x = x(\[Gamma]) : \[Gamma] (phase) must be within a particular range so that the other entries remain unimodular *)
% (* x = -0.9481448970234081` - 0.31783840902016647`*I; first try for phase = 0.551478475784210 (2018-01-01) *)
% x1 = Exp[I Pi Random[Real, {0.38629047387009, 0.8301134750001507}]];
% x2 = Exp[I Pi Random[Real, {0.9700212764097545, 1.413709525400110}]];
% t = Random[Integer];
% (* take either of two branches... *)
% x = t*x1 + (1 - t)*x2
% 
% 
% (* get T8 = T8(x) with y = y(x), z = z(x), u = u(x) *)
% Subscript[\[Chi], 1] = 5 - 3 Sqrt[5] + 4 Sqrt[5 - 2 Sqrt[5]] + (1 + I) 1/x + Subscript[\[Zeta], 2];
% Subscript[\[Chi], 2] = 45 - 9 I - (17 - 7 I) Sqrt[5] - 2 Sqrt[620 - 500 I - (285 - 234 I) Sqrt[5]] - (1 + I) Subscript[\[Chi], 1]/x;
% Subscript[\[Chi], 3] = 10 I - 10 + (6 - 2 I) Sqrt[5] + 4 Sqrt[10 - 10 I - (5 - 4 I) Sqrt[5]] + (1 + I) Subscript[\[Chi], 2]/x;
% Subscript[\[Chi], 4] = Sqrt[5] - 1 + I Sqrt[2] Sqrt[5 + Sqrt[5]] + (1 + I) Subscript[\[Chi], 3]/x;
% Subscript[\[Chi], 5] = (1 + I) Sqrt[1 - Sqrt[5] + Subscript[\[Zeta], 4]] - (1 - I) (2 x^2 + Conjugate[Sqrt[Subscript[\[Chi], 4]]]) - x (1 - 4 I + Sqrt[5] + Subscript[\[Zeta], 1]);
% Subscript[\[Chi], 6] = Sqrt[5] - 1 - I - Sqrt[5 - 2 Sqrt[5]] + (1 - I) x;
% y = 1/4 Subscript[\[Chi], 5]/Subscript[\[Chi], 6];
% z = x/(x - y) ((1 + I) (-1)^(4/5) - x - y - 1 + I);
% u = I/z (x^2 - z x + x y + y z)/(1 + I - (-1)^(3/10) - (-1)^(4/5));
% T = ( {
%     {0, 0, 0, 0, 0, 0, 0, 0},
%     {0, 0, 0, 0, 0, 3, 15, 18},
%     {0, 8, 8, 8, 8, 13, 15, 8},
%     {0, 2, 2, 2, 2, 17, 5, 2},
%     {0, 0, 0, 0, 0, 7, 5, 12},
%     {0, 0, 0, 0, 0, 0, 10, 10},
%     {0, 0, 0, 0, 0, 10, 0, 10},
%     {0, 0, 0, 0, 0, 10, 10, 0}
%    } );
% T8 = Exp[I Pi T 1/10]*( {
%      {1, 1, 1, 1, 1, 1, 1, 1},
%      {1, x, y, z, -((y z)/x), 1, 1, 1},
%      {1, -(u/y), u/x, -((u z)/(x y)), -((u z)/x^2), 1, 1, 1},
%      {1, y, x , -((y z)/x), z , 1, 1, 1},
%      {1, u/x, -(u/y), -((u z)/x^2), -((u z)/(x y)), 1, 1, 1},
%      {1, -I x/z, I x/z, -I z/x, I z/x, 1, 1, 1},
%      {1, u, -u, -((u z^2)/x^2), (u z^2)/x^2, 1, 1, 1},
%      {1, I (u z)/x, I (u z)/x, -I (u z)/x, -I (u z)/x, 1, 1, 1}
%     } );
% 
% 
% (* due to the construction Arg[u]/2./Pi should be in [-0.5, 0.5] *)
% Arg[u]/2./Pi
% T8 // MatrixForm // N // Chop
% 
% 
% (* check CHM class *)
% T8.ConjugateTranspose[T8] // N // MatrixForm // Chop
% Abs[T8] // MatrixForm // N // Chop
% 
% 
% (* [x] = 5/10 or 10/10 or 13/10, defect = 7, BH(8, 20) *)
% (* three equivalent Butsons (each has different order of columns {2,3,4,5}) *)
% B8A = Exp[I Pi 1/10 ( {
%       {0, 0, 0, 0, 0, 0, 0, 0},
%       {0, 5, 10, 13, 8, 3, 15, 18},
%       {0, 10, 5, 18, 3, 13, 15, 8},
%       {0, 12, 7, 10, 15, 17, 5, 2},
%       {0, 17, 2, 15, 10, 7, 5, 12},
%       {0, 7, 17, 3, 13, 0, 10, 10},
%       {0, 2, 12, 8, 18, 10, 0, 10},
%       {0, 15, 15, 5, 5, 10, 10, 0}
%      } )];
% 
% 
% (* [x] = 8/10, defect = 11, BH(8, 20) *)
% B8B = Exp[I Pi 1/10 ( {
%       {0, 0, 0, 0, 0, 0, 0, 0},
%       {0, 8, 10, 13, 5, 3, 15, 18},
%       {0, 18, 10, 3, 5, 13, 15, 8},
%       {0, 12, 10, 7, 15, 17, 5, 2},
%       {0, 2, 10, 17, 15, 7, 5, 12},
%       {0, 10, 0, 0, 10, 0, 10, 10},
%       {0, 10, 0, 10, 0, 10, 0, 10},
%       {0, 0, 0, 10, 10, 10, 10, 0}
%      } )];
% ------------------------------------------------------------------------------
% :
% Module[{f = Interpolation[V, InterpolationOrder -> 20][x], g, h},
%  g = D[f, {x, 5}];
%  h = Convolve[g, Exp[-100 x^2], x, y];
%  Print[h];
% Plot[{g}, {x, 0.02, 0.45}]]
% Convolve[InterpolatingFunction[{{0.01, 0.479}}, <>][x], E^(-100 x^2), x, y]
% 
% VV = If[#2 > .5, List@##, {#1, #2 + 1, #3}] & @@@ UU;
% WW = ({#.norm3, #.norm2} &) /@ VV;
% ------------------------------------------------------------------------------
