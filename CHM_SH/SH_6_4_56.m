function Y = SH_6_4_56
% ------------------------------------------------------------------------------
% 2024-08-30
%
% symmetric CHM(6) with d = 4 and L = 56
% ------------------------------------------------------------------------------

    a = exp(1j*(atan(2*sqrt(2)-sqrt(5))-pi));
    b = exp(1j*(pi-acot(sqrt(2))));
    c = exp(1j*atan(2*sqrt(2)+sqrt(5)));
    d = exp(1j*atan((sqrt(3)-2*sqrt(2))/(1+2*sqrt(6))));
    e = -b;
    x = exp(-2j*pi/3);
    y = 1/x;
    f = x*d;

    Y = [
        1  1  1  1  1  1;
        1  x  a  b  c  d;
        1  a  y  e  f  c;
        1  b  e -1  e  b;
        1  c  f  e  y  a;
        1  d  c  b  a  x;
    ];



end
