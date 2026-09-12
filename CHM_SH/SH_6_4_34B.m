function Y = SH_6_4_34B
% ------------------------------------------------------------------------------
% 2024-08-30
% variant of SH_6_4_34.m
% ------------------------------------------------------------------------------

%    b=0.3468675145425457 -     0.9379141364512409j; % original value
%    but Math. cannot solve the associated equation -- it finds a different value:

    b=-(1/2)+1/sqrt(2)-1/2*1j*sqrt(1+2*sqrt(2));
    a=(-b*(1+b)^2+sqrt(2)*sqrt(b^2*(1+2*b+2*b^3+b^4)))/(-1+b*(2+b));
    c=(-a-a^2-2*a*b)/(a+b^2);
    x=-1-a-2*b-c;
    d=a/b;
    f=a/c;
    y=a/x;
    Y=[1,  1,  1,  1,  1,  1;
       1,  a,  b,  b,  c,  x;
       1,  d, -1, -d, -c,  c;
       1,  d, -b, -1, -d,  b;
       1,  f, -f, -b, -1,  b;
       1,  y,  f,  d,  d,  a];

end
