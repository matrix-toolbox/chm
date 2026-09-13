function C8R_check(z)
% 2021-06-24
% https://homepages.math.uic.edu/~jan/FactorBench/cyclic8e1.html (accessed: 2021-06-24)
%
%
% Y = sinkhorn(C8_CSEED) ; circulant seed
% Y = dephase(Y)
% z = Y(1, :);
% C8R_check(z)


% some C-8-R solution:
%
% z(1) =  5.18351470931885E-01  +3.14750403952489E-01*1j;
% z(2) =  7.38451529753908E-01  -8.25057190030791E-01*1j;
% z(3) = -3.09651833066213E-02  -6.44561566127122E-01*1j;
% z(4) =  7.76221780791704E-01  -2.17337913503343E+00*1j;
% z(5) = -7.38451529753908E-01  +8.25057190030791E-01*1j;
% z(6) = -5.18351470931885E-01  -3.14750403952490E-01*1j;
% z(7) = -7.76221780791704E-01  +2.17337913503342E+00*1j;
% z(8) =  3.09651833066214E-02  +6.44561566127122E-01*1j;

z0 = z(1);
z1 = z(2);
z2 = z(3);
z3 = z(4);
z4 = z(5);
z5 = z(6);
z6 = z(7);
z7 = z(8);


% Stupid solution:
% Actually this is not stupid, but just wrong.
% The following formulas are taken from the pattern of V_N
% just to confirm that this has nothing to do with cyclic-N-roots...
% update @2023-02-07
% Really?!

s1 = z0/z7 + z1/z0 + z2/z1 + z3/z2 + z4/z3 + z5/z4 + z6/z5 + z7/z6;
s2 = z0/z6 + z1/z7 + z2/z0 + z3/z1 + z4/z2 + z5/z3 + z6/z4 + z7/z5;
s3 = z0/z5 + z1/z6 + z2/z7 + z3/z0 + z4/z1 + z5/z2 + z6/z3 + z7/z4;
s4 = z0/z4 + z1/z5 + z2/z6 + z3/z7 + z4/z0 + z5/z1 + z6/z2 + z7/z3;
s5 = z0/z3 + z1/z4 + z2/z5 + z3/z6 + z4/z7 + z5/z0 + z6/z1 + z7/z2;
s6 = z0/z2 + z1/z3 + z2/z4 + z3/z5 + z4/z6 + z5/z7 + z6/z0 + z7/z1;
s7 = z0/z1 + z1/z2 + z2/z3 + z3/z4 + z4/z5 + z5/z6 + z6/z7 + z7/z0;



OUR_EQUATION=[s1;s2;s3;s4;s5;s6;s7]


# Correct solution:
r1 = z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7;
r2 = z0*z1 + z1*z2 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7 + z7*z0;
r3 = z0*z1*z2 + z1*z2*z3 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7 + z6*z7*z0 + z7*z0*z1;
r4 = z0*z1*z2*z3 + z1*z2*z3*z4 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7 + z5*z6*z7*z0 + z6*z7*z0*z1 + z7*z0*z1*z2;
r5 = z0*z1*z2*z3*z4 + z1*z2*z3*z4*z5 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7 + z4*z5*z6*z7*z0 + z5*z6*z7*z0*z1 + z6*z7*z0*z1*z2 + z7*z0*z1*z2*z3;
r6 = z0*z1*z2*z3*z4*z5 + z1*z2*z3*z4*z5*z6 + z2*z3*z4*z5*z6*z7 + z3*z4*z5*z6*z7*z0 + z4*z5*z6*z7*z0*z1 + z5*z6*z7*z0*z1*z2 + z6*z7*z0*z1*z2*z3 + z7*z0*z1*z2*z3*z4;
r7 = z0*z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z6*z7 + z2*z3*z4*z5*z6*z7*z0 + z3*z4*z5*z6*z7*z0*z1 + z4*z5*z6*z7*z0*z1*z2 + z5*z6*z7*z0*z1*z2*z3 + z6*z7*z0*z1*z2*z3*z4 + z7*z0*z1*z2*z3*z4*z5;
r8 = z0*z1*z2*z3*z4*z5*z6*z7-1

CYCLIC_N_ROOTS=[r1;r2;r3;r4;r5;r6;r7;r8]
