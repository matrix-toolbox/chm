function V = LN15cos()
% 2023-01-11
%
% Static example for k = 3 => N = 15.

    clc
    k = 3;
    N = 3 + 4*k;
    printf('Wait, calculating solution...');
    do
        [a info]=fsolve(@uc, 2*pi*rand(1, 2*k+1));
        for j=1:128
            [a info]=fsolve(@uc, a);
        end
        s = sum(info);
        printf('.');
    until s < 1e-13

    printf(' Solved!\n');
    L = [
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 a(1) -a(1) a(2) -a(2) a(3) -a(3) a(4) -a(4) a(5) -a(5) a(6) -a(6) a(7) -a(7);
        0 -a(1) a(1) -a(2) a(2) -a(3) a(3) -a(4) a(4) -a(5) a(5) -a(6) a(6) -a(7) a(7);
        0 a(7) -a(7) a(1) -a(1) a(2) -a(2) a(3) -a(3) a(4) -a(4) a(5) -a(5) a(6) -a(6);
        0 -a(7) a(7) -a(1) a(1) -a(2) a(2) -a(3) a(3) -a(4) a(4) -a(5) a(5) -a(6) a(6);
        0 a(6) -a(6) a(7) -a(7) a(1) -a(1) a(2) -a(2) a(3) -a(3) a(4) -a(4) a(5) -a(5);
        0 -a(6) a(6) -a(7) a(7) -a(1) a(1) -a(2) a(2) -a(3) a(3) -a(4) a(4) -a(5) a(5);
        0 a(5) -a(5) a(6) -a(6) a(7) -a(7) a(1) -a(1) a(2) -a(2) a(3) -a(3) a(4) -a(4);
        0 -a(5) a(5) -a(6) a(6) -a(7) a(7) -a(1) a(1) -a(2) a(2) -a(3) a(3) -a(4) a(4);
        0 a(4) -a(4) a(5) -a(5) a(6) -a(6) a(7) -a(7) a(1) -a(1) a(2) -a(2) a(3) -a(3);
        0 -a(4) a(4) -a(5) a(5) -a(6) a(6) -a(7) a(7) -a(1) a(1) -a(2) a(2) -a(3) a(3);
        0 a(3) -a(3) a(4) -a(4) a(5) -a(5) a(6) -a(6) a(7) -a(7) a(1) -a(1) a(2) -a(2);
        0 -a(3) a(3) -a(4) a(4) -a(5) a(5) -a(6) a(6) -a(7) a(7) -a(1) a(1) -a(2) a(2);
        0 a(2) -a(2) a(3) -a(3) a(4) -a(4) a(5) -a(5) a(6) -a(6) a(7) -a(7) a(1) -a(1);
        0 -a(2) a(2) -a(3) a(3) -a(4) a(4) -a(5) a(5) -a(6) a(6) -a(7) a(7) -a(1) a(1);
    ];
    K = exp(1j*L);
    V = dephase(K);
    isBH(V, 10000);
    printf('d(V) = %d\n', ud(V, 'S', 1e-8));
    printf('#Lambda(V) = ... (not calculated, please do it manually)\n');
    A = mod(angle(V) + 2 * pi, 2 * pi) / 2 / pi;
    save('LN_cos.data', 'A');

end

function y=uc(a)
 y(1)= 1/2 + cos(a(1))+cos(a(2))+cos(a(3))+cos(a(4))+cos(a(5))+cos(a(6))+cos(a(7));
 y(2)= 1/2 + cos(2*a(1))+cos(2*a(2))+cos(2*a(3))+cos(2*a(4))+cos(2*a(5))+cos(2*a(6))+cos(2*a(7));
 y(3) = 1/2 + cos(a(1) + a(2)) + cos(a(2) + a(3)) + cos(a(3) + a(4)) + cos(a(4) + a(5)) + cos(a(5) + a(6)) + cos(a(6) + a(7)) + cos(a(7) + a(1));
 y(4) = 1/2 + cos(a(1) + a(3)) + cos(a(2) + a(4)) + cos(a(3) + a(5)) + cos(a(4) + a(6)) + cos(a(5) + a(7)) + cos(a(6) + a(1)) + cos(a(7) + a(2));
 y(5) = 1/2 + cos(a(1) + a(4)) + cos(a(2) + a(5)) + cos(a(3) + a(6)) + cos(a(4) + a(7)) + cos(a(5) + a(1)) + cos(a(6) + a(2)) + cos(a(7) + a(3));
 y(6) = 1/2 + cos(a(1) - a(2)) + cos(a(2) - a(3)) + cos(a(3) - a(4)) + cos(a(4) - a(5)) + cos(a(5) - a(6)) + cos(a(6) - a(7)) + cos(a(7) - a(1));
 y(7) = 1/2 + cos(a(1) - a(3)) + cos(a(2) - a(4)) + cos(a(3) - a(5)) + cos(a(4) - a(6)) + cos(a(5) - a(7)) + cos(a(6) - a(1)) + cos(a(7) - a(2));
 y(8) = 1/2 + cos(a(1) - a(4)) + cos(a(2) - a(5)) + cos(a(3) - a(6)) + cos(a(4) - a(7)) + cos(a(5) - a(1)) + cos(a(6) - a(2)) + cos(a(7) - a(3));
end
