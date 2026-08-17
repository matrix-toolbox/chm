function Y = SG_9_0_18
% ------------------------------------------------------------------------------
% 2022-12-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Self G-dual CHM with very nice parameters.
% Looks like N9_0...
% ------------------------------------------------------------------------------

    y = -1/4 + 1j*sqrt(15)/4;
    w = exp(2j*pi/3);
    z = 7/128 + 1j*sqrt((128^2-49)/128^2); % beautiful number :)

    N9_0 = [
        1    1       1       1     1     1      1     1     1    ;
        1    y       y^2    -1    -y     y      y^3   y^3   y    ;
        1    y^2     y^4    -y     y^2  -y^3    y^4   y^2   1    ;
        1   -y^4    -y^3     y^3  -y^3  -y^2   -y^4  -y^2  -1    ;
        1   -1       y^2    -1/y   1     1      y     y^2   1/y  ;
        1    y^3    -y      -1     y^2   y^3    y     y     y    ;
        1    y       1      -1/y   y^2   y^2   -1     1     1/y  ;
        1    y^4     y^2    -y     y^4   y^2    y^2  -y^3   1    ;
        1    y^3     y^4    -y^3   y^4   y^2    y^3   y^2  -1    ;
    ];

    DL = diag([1 1 1 1 -y^4 -y^3 1 y 1]);
    DR = diag([1 1 1 1 -1 -y -y^3, z, z*y]);

    Y = DL*N9_0*DR;

end




