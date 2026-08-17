function [iiMax, nP, ZTolerance, YPattern, muFactor] = SH_7_0_97
% ------------------------------------------------------------------------------
% 2023-02-02 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Symmetric CHM of size N = 7 with generic d = 0.
% Convergence only for particular starting points...
% This pattern contains solutions with:
% 1) #L = 6  ==> P7 in BH(7, 6) with d = 3
% 2) #L = 7  ==> F7 in BH(7, 7)
% 3) #L = 97 ==> Q7
% ------------------------------------------------------------------------------
%     H = exp(2j*pi*[
%        0 0 0 0 0 0 0;
%        0 1 5 3 5 2 3;
%        0 5 1 5 3 3 2;
%        0 3 5 1 2 5 3;
%        0 5 3 2 1 3 5;
%        0 2 3 5 3 1 5;
%        0 3 2 3 5 5 1;
%    ]/6);
%
%    H = exp(2j*pi*[
%        0 0 0 0 0 0 0;
%        0 4 2 5 5 2 2;
%        0 2 4 2 5 5 2;
%        0 5 2 2 3 0 4;
%        0 5 5 3 2 3 1;
%        0 2 5 0 3 2 4;
%        0 2 2 4 1 4 5;
%    ]/6);
% ------------------------------------------------------------------------------

    iiMax = 40000;
    nP = 6;
    ZTolerance = 5e-13;
    YPattern = @pattern;
    muFactor = 0.001;
end

function Y = pattern(p)
    Y = [
        1  1     1     1     1     1     1   ;
        1  p(1)  p(2)  p(3)  p(4)  p(5)  p(6);
        1  p(2)  p(3)  p(5)  p(6)  p(4)  p(1);
        1  p(3)  p(5)  p(4)  p(1)  p(6)  p(2);
        1  p(4)  p(6)  p(1)  p(3)  p(2)  p(5);
        1  p(5)  p(4)  p(6)  p(2)  p(1)  p(3);
        1  p(6)  p(1)  p(2)  p(5)  p(3)  p(4);
    ];
end
