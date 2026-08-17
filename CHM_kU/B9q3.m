function B = B9q3(j)
% ------------------------------------------------------------------------------
% 2022-12-25 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Taken from Butson Home;
% https://wiki.aalto.fi/display/Butson/Matrices+up+to+monomial+equivalence
% ------------------------------------------------------------------------------

    BH{1} = [ 0 0 0 0 0 0 0 0 0;0 0 0 1 1 1 2 2 2;0 0 0 2 2 2 1 1 1;0 1 2 0 1 2 0 1 2;0 1 2 1 2 0 2 0 1;0 1 2 2 0 1 1 2 0;0 2 1 0 2 1 0 2 1;0 2 1 1 0 2 2 1 0;0 2 1 2 1 0 1 0 2];
    BH{2} = [ 0 0 0 0 0 0 0 0 0;0 0 0 1 1 1 2 2 2;0 0 0 2 2 2 1 1 1;0 1 2 0 1 2 0 1 2;0 1 2 1 2 0 2 0 1;0 1 2 2 0 1 1 2 0;0 2 1 0 2 1 1 0 2;0 2 1 1 0 2 0 2 1;0 2 1 2 1 0 2 1 0];
    BH{3} = [ 0 0 0 0 0 0 0 0 0;0 0 0 1 1 1 2 2 2;0 0 0 2 2 2 1 1 1;0 1 2 0 1 2 0 1 2;0 1 2 1 2 0 2 0 1;0 1 2 2 0 1 1 2 0;0 2 1 0 2 1 2 1 0;0 2 1 1 0 2 1 0 2;0 2 1 2 1 0 0 2 1];

    MAX_BH = size(BH, 2);
    if j > MAX_BH
        error("Index out of scope! Must be in range 1 .. %d", MAX_BH);
    else
        B = exp(2j * pi * BH{j} / 3);
    end

end

