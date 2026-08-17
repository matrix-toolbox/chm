function H = HH_8_7_18(p)
% ------------------------------------------------------------------------------
% 2023-01-14 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% A 2-parametric Hermitian family of CHM of order N = 8 stemming from BH(8, 2).
% Generic defect = 7, #L = 18.
%
% p1   p2     d    #L     q : BH(8, q)
% ====================================
% 4/5  3/5    7    10     10          
% 1/4  2/3    7    12     12          
% r1   r2     7*   18*     -           <-- generic values
% 1/4  1/4    9     4      4               for random pj = rj in [0, 1)
% 1/3  1/3    9     6      6          
% 3/5  3/5    9    10     10          
% 1/2  3/4    11    4      4          
% r    0      11    6      -          
% 0    0      21    2      2          
% ------------------------------------------------------------------------------

    b = exp(2j*pi * p(1));
    c = exp(2j*pi * p(2));

    H = [
        1   1   1   1   1   1   1   1 ;
        1  -1   b   c  -1   1  -c  -b ;
        1   b' -1  -b'  1   b' -b' -1 ;
        1   c' -b   1  -1  -c' -1   b ;
        1  -1   1  -1   1  -1  -1   1 ;
        1   1   b  -c  -1  -1   c  -b ;
        1  -c' -b  -1  -1   c'  1   b ;
        1  -b' -1   b'  1  -b'  b' -1 ;
    ];

end
