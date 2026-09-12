function H = M6(p)
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric (by construction) CHM found by M. Matolcsi and F. Szollosi (2006);
%
% >> H = M6(rand);
%
% Generic values of parameter p1 provide d=4 and #L=58.
% M6(1) in BH(6, 6).
% ------------------------------------------------------------------------------

    p = pi*(p + 1)/2; % single parameter in (pi/2, pi] or (3*pi/2, 2*pi]
    if rand < 0.5
        p = pi + p;
    end

    a = exp(2j * pi * p);

    t1 = a^2 + 2 * a - 1;
    t2 = 1 + a^2;
    t3 = 1 - 6 * a^2 + a^4;
    t4 = sqrt(1 + 4 * a^2 - 58 * a^4 + 4 * a^6 + a^8);
    t5 = sqrt(16 - abs(t2)^2) / (4 * abs(t2));
    t6 = 1j * t1 * sqrt(16 - abs(t1)^2) / (4 * abs(t1));

    b = (t3 - t4) / (4 * t1);
    c = (t3 + t4) / (4 * t1);
    d = t2 / 4 + 1j * t2 * t5;
    e = t2 / 4 - 1j * t2 * t5;
    f = t1 / 4 + t6;
    g = t1 / 4 - t6;

    H = [
        1,  1,  1,  1,  1,  1;
        1, -1,  a,  a, -a, -a;
        1,  a,  b,  c, -d, -e;
        1,  a,  c,  b, -e, -d;
        1, -a, -d, -e,  f,  g;
        1, -a, -e, -d,  g,  f;
    ];

end
