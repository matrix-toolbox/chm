function H = S8(p)
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 4-parametric family of CHM stemming from BH(8, 4) found by M. Matolcsi et al.
% https://arxiv.org/abs/quant-ph/0607073
%
% Matolcsi et al. [arXiv:0607073] observe that S8 and F8 intersect at H8 for some parameters.
% S8 and F8 both have defect 5 but they are not equivalent,
% hence they are not subfamilies of a common family with defect > 5.
%
% >> S8(rand(1, 4));
% ------------------------------------------------------------------------------

    S8 = exp(2j*pi*[
        0 0 0 0 0 0 0 0;
        0 0 2 2 2 1 3 0;
        0 1 1 3 0 2 2 3;
        0 1 3 1 2 3 1 3;
        0 2 3 1 0 1 3 2;
        0 2 1 3 2 0 0 2;
        0 3 2 2 0 3 1 1;
        0 3 0 0 2 2 2 1;
    ]/4);

    R = [
        0 0         0         0         0 0         0         0        ;
        0 p(4)      p(4)      p(4)      0 p(3)+p(4) p(3)+p(4) p(4)     ;
        0 p(1)-p(4) p(2)-p(4) p(2)-p(4) 0 0         0         p(1)-p(4);
        0 p(1)      p(2)      p(2)      0 p(3)+p(4) p(3)+p(4) p(1)     ;
        0 0         p(2)-p(4) p(2)-p(4) 0 p(3)      p(3)      0        ;
        0 p(4)      p(2)      p(2)      0 p(4)      p(4)      p(4)     ;
        0 p(1)-p(4) 0         0         0 p(3)      p(3)      p(1)-p(4);
        0 p(1)      p(4)      p(4)      0 p(4)      p(4)      p(1)     ;
    ];

    H = S8 .* exp(2j * pi * R);

end
