function H = BH10_6(p)
% ------------------------------------------------------------------------------
% 2008-05-26 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% An example of 1-parametric affine family stemming from BH(10, 6).
% More Butson-type matrices can be found in the "Butson Home" at:
% https://wiki.aalto.fi/display/Butson/Butson+Home
%
% >> H = BH10_6(rand);
% ------------------------------------------------------------------------------

    Z10 = exp(2j*pi*[
        0 0 0 0 0 0 0 0 0 0;
        0 3 1 1 4 3 4 0 1 4;
        0 3 0 4 5 2 0 2 4 2;
        0 0 5 3 3 3 1 3 1 5;
        0 1 3 5 0 1 4 3 2 4;
        0 3 3 1 4 0 2 3 5 0;
        0 3 3 4 2 5 0 5 1 2;
        0 1 4 2 0 4 2 0 4 3;
        0 5 0 1 2 5 4 3 3 2;
        0 5 2 4 2 2 3 0 4 0;
    ]/6);

    R = [
        0 0 0 0 0 0 0 0 0 0;
        0 0 p p p 0 p 0 p p;
        0 0 0 p 0 0 p 0 0 p;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 p 0 0 0 p p;
        0 0 p p p 0 p 0 p p;
        0 0 0 p 0 0 p 0 0 p;
        0 0 p 0 p 0 p 0 0 0;
        0 0 p 0 p 0 p 0 0 0;
        0 0 0 0 p 0 0 0 p p;
    ];

    H = Z10 .* exp(2j * pi * R);

end
