function H = A8
% ------------------------------------------------------------------------------
% 2016-05-31 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric isolated CHM of order N = 8 found by WB;
% [1] https://link.springer.com/article/10.1007/s11786-018-0379-8
% [2] https://arxiv.org/pdf/2204.11727.pdf
%
% >> H = A8;
% ------------------------------------------------------------------------------

    a = 1/3 + 1j*sqrt(8)/3;

    H = [
        1    1    1    1    1    1    1    1   ;
        1 -a^2   -1    a^2  a   -a   -a    a   ;
        1   -1   -1   -a    a'   1   -a'   a   ;
        1  a^2   -a    a   -1    a   -a   -a^2 ;
        1    a    a'  -1    1   -1   -a'  -a   ;
        1   -a    1    a   -1    a'  -a'  -1   ;
        1   -a   -a'  -a   -a'  -a'   1   -a   ;
        1    a    a   -a^2 -a   -1   -a    a^2 ;
    ];

end
