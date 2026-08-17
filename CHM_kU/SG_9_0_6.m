function Y = SG_9_0_6
% ------------------------------------------------------------------------------
% 2022-12-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% (First) self G-dual matrix B from BH(9, 6) with vanishing defect d(B) = 0.
% ------------------------------------------------------------------------------

    Y = exp(2j*pi*[
        0   0   0   0   0   0   0   0   0;
        0   5   3   2   5   3   2   1   5;
        0   3   3   0   1   5   4   1   3;
        0   2   0   2   0   2   4   4   4;
        0   5   1   0   3   3   4   3   1;
        0   3   5   2   3   5   2   5   1;
        0   0   2   4   2   0   2   4   4;
        0   3   3   4   5   1   0   3   1;
        0   1   5   4   3   3   0   1   3;
    ]/6);

    w = exp(2j*pi/3);

    DL=diag([1,1,1,1,w,1,1,1,w]);
    DR=diag([1,1,1,1,w,1,1,w,w^2]);

    Y = DL * Y * DR;
end

