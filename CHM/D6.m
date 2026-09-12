function H = D6(p)
% ------------------------------------------------------------------------------
% 2006-12-04 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Example of an affine family depending on 1 parameter, stemming from BH(6, 4);
%
% >> H = D6(rand);
% ------------------------------------------------------------------------------

    D6 = exp(2j*pi*[
        0 0 0 0 0 0;
        0 2 1 3 3 1;
        0 1 2 1 3 3;
        0 3 1 2 1 3;
        0 3 3 1 2 1;
        0 1 3 3 1 2;
    ]/4);

    R = [
        0  0  0  0  0  0;
        0  0  0  0  0  0;
        0  0  0  p  p  0;
        0  0 -p  0  0 -p;
        0  0 -p  0  0 -p;
        0  0  0  p  p  0;
    ];

    H = D6 .* exp(2j * pi * R);

end
