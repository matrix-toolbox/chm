function B = B9q10
% ------------------------------------------------------------------------------
% 2022-09-19 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Taken from Butson Home;
% https://wiki.aalto.fi/display/Butson/Matrices+up+to+monomial+equivalence
% ------------------------------------------------------------------------------

    LB = [
         0 0 0 0 0 0 0 0 0;
         0 0 1 2 4 5 6 6 8;
         0 1 4 7 9 3 5 9 5;
         0 2 7 3 8 8 2 4 6;
         0 4 9 8 4 6 2 8 3;
         0 5 3 8 6 2 0 4 8;
         0 6 5 2 2 0 8 7 4;
         0 6 9 4 8 4 7 2 2;
         0 8 5 6 3 8 4 2 0];

    B = exp(2j * pi * LB / 10); % in BH(9, 10) there is only one element

end

