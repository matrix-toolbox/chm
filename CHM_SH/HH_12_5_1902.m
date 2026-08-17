function [iiMax, nP, ZTolerance, YPattern, muFactor] = HH_12_5_1902()
% ------------------------------------------------------------------------------
% 2023-01-23 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian CHM of size N = 12 with generic d = 5 and #L = 1902.
% Very slow convergence only for particular seeds.
% ------------------------------------------------------------------------------
% Generic values: d = 5, #L = 1902 were observed only ONCE! :)
% See: "HH_12_5_1902_20230123T121254.data".
%
% This pattern contains many other solutions:
% d = 9, #L = 272 ---> see HH_12_9_272.m
% d = 9, #L = 294 ---> see HH_12_9_294.m
% d = 9, #L = 314 ---> see HH_12_9_314.m
% d = 9, #L = 328 ---> see HH_12_9_328.m
% d = 9, #L = 622 ---> see HH_12_9_622.m
% d =11, #L = 6   ---> see HH_12_11_6.m  : BH(12, 6)
% d =13, #L = 18  ---> see HH_12_13_18.m : 2-parametric family
% ------------------------------------------------------------------------------

    iiMax = 50000;
    nP = 50;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.0004;
end

function Y = pattern(p)
    Y = [
        1   1       1       1       1       1       1       1      1        1       1       1      ;
        1   1       p(1)    p(2)    p(3)    p(4)    p(5)    p(6)    p(7)    p(8)    p(9)    p(10)  ;
        1   p(1)'   1       p(11)'  p(12)'  p(13)'  p(14)'  p(15)' -p(16)   p(17)'  p(18)'  p(19)' ;
        1   p(2)'   p(11)   1       p(20)   p(21)   p(22)'  p(23)'  p(24)'  p(25)'  p(26)   p(27)  ;
        1   p(3)'   p(12)   p(20)'  1       p(28)   p(29)'  p(30)'  p(31)'  p(32)'  p(33)   p(34)  ;
        1   p(4)'   p(13)   p(21)'  p(28)'  1       p(35)'  p(36)' -p(37)'  p(38)'  p(39)   p(40)  ;
        1   p(5)'   p(14)   p(22)   p(29)   p(35)  -1       p(41)  -p(7)    p(42)   p(43)   p(44)  ;
        1   p(6)'   p(15)   p(23)   p(30)   p(36)   p(41)' -1       p(16)   p(45)   p(46)   p(47)  ;
        1   p(7)'  -p(16)'  p(24)   p(31)  -p(37)  -p(7)'   p(16)' -1       p(37)  -p(24)  -p(31)  ;
        1   p(8)'   p(17)   p(25)   p(32)   p(38)   p(42)'  p(45)'  p(37)' -1       p(48)   p(49)  ;
        1   p(9)'   p(18)   p(26)'  p(33)'  p(39)'  p(43)'  p(46)' -p(24)'  p(48)' -1       p(50)  ;
        1   p(10)'  p(19)   p(27)'  p(34)'  p(40)'  p(44)'  p(47)' -p(31)'  p(49)'  p(50)' -1      ;
    ];
end
