function Y = SH_6_4_48_DSC
% ------------------------------------------------------------------------------
% 2024-08-26
%
% symmetric CHM(6) with d = 4 and L = 48 and doubly symmetric core
% ------------------------------------------------------------------------------

    c = exp(-2j*atan(sqrt(7+2*sqrt(13))));
    d = exp(-1j*atan(sqrt(18/(1+5*sqrt(13)-8*sqrt(2+sqrt(13))))));
    e = exp(1j*(pi-atan(sqrt(23+7*sqrt(13))/2/sqrt(2))));

    b = c^2/d;
    g = -c^2;
    f = -c^2/e;

    Y = [
        1  1  1  1  1  1;
        1  1  b  c  d  e;
        1  b  f -c  g  d;
        1  c -c -1 -c  c;
        1  d  g -c  f  b;
        1  e  d  c  b  1;
    ];



end
