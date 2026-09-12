function H = H36L
% ------------------------------------------------------------------------------
% 2025-04-14
% https://webspace.maths.qmul.ac.uk/l.h.soicher/designtheory.org/library/encyc/topics/had.pdf
% ------------------------------------------------------------------------------

% Define a Latin square of order 6
L = [1 2 3 4 5 6;
     2 3 4 5 6 1;
     3 4 5 6 1 2;
     4 5 6 1 2 3;
     5 6 1 2 3 4;
     6 1 2 3 4 5];

% List of cells in the 6x6 array
cells = zeros(36, 2);
index = 1;
for i = 1:6
    for j = 1:6
        cells(index, :) = [i, j];
        index = index + 1;
    end
end

% Initialize Hadamard matrix H
H = -ones(36, 36);  % Fill with -1s

% Fill H according to the rule
for a = 1:36
    i1 = cells(a, 1);
    j1 = cells(a, 2);
    v1 = L(i1, j1);
    for b = 1:36
        if a == b
            continue;  % Keep diagonal as -1
        end
        i2 = cells(b, 1);
        j2 = cells(b, 2);
        v2 = L(i2, j2);
        if (i1 == i2) || (j1 == j2) || (v1 == v2)
            H(a, b) = 1;
        end
    end
end

    printf("Is Hadamard matrix? %d\n", isequal(H * transpose(H), 36 * eye(36)));


end

