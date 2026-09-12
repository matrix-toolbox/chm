function H = P7(p)
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 2-parametric CHM of order N = 7 stemming from BH(7, 6), found by M. Petrescu;
%
% >> H = P7(rand);
% ------------------------------------------------------------------------------

    P7 = exp(2j*pi*[
        0 0 0 0 0 0 0;
        0 1 4 5 3 3 1;
        0 4 1 3 5 3 1;
        0 5 3 1 4 1 3;
        0 3 5 4 1 1 3;
        0 3 3 1 1 4 5;
        0 1 1 3 3 5 4;
    ]/6);

    R = [
        0 0 0  0  0 0 0;
        0 p p  0  0 0 0;
        0 p p  0  0 0 0;
        0 0 0 -p -p 0 0;
        0 0 0 -p -p 0 0;
        0 0 0  0  0 0 0;
        0 0 0  0  0 0 0;
    ];

    H = P7 .* exp(2j * pi * R);

end
