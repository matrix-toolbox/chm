function Y = BH_6_4_6;
% ------------------------------------------------------------------------------
% 2023-02-02 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric BH(6, 6) with d = 4 and #L = 6.
% ------------------------------------------------------------------------------

    a = exp(2j*pi/6);
    Y = [
        1   1   1   1   1   1  ;
        1   a  -1  -a'  a' -a  ;
        1  -1  -1   1  -1   1  ;
        1  -a'  1  -a  -a  -a' ;
        1   a' -1  -a   a  -a' ;
        1  -a   1  -a' -a' -a  ;
    ];

end
