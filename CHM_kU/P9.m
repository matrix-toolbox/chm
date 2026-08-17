function P = P9
% ------------------------------------------------------------------------------
% 2022-09-29 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Matrix representation of AME(4, 3) state as a permutation matrix P in P(9).
% ------------------------------------------------------------------------------

    I3 = eye(3);
    I9 = kron(I3, I3);
    w = exp(2j * pi / 3);
    F3 = [
        1 1 1;
        1 w w^2;
        1 w^2 w
    ]/sqrt(3);

    ket0 = [1; 0; 0];
    ket1 = [0; 1; 0];
    ket2 = [0; 0; 1];
    p0 = ket0 * ket0';
    p1 = ket1 * ket1';
    p2 = ket2 * ket2';

    X1=[
        0 0 1;
        1 0 0;
        0 1 0];
    X2=[
        0 1 0;
        0 0 1;
        1 0 0;
    ];


    A = ...
    % kron(I3,I3,F3,F3) *...
    %(kron(I3,I3,I3,p0) + kron(I3,X1,I3,p1) + kron(I3,X2,I3,p2)) * ...
    %(kron(I3,I3,I3,p0) + kron(X1,I3,I3,p1) + kron(X2,I3,I3,p2)) * ...
    %(kron(I3,I3,p0,I3) + kron(I3,X1,p1,I3) + kron(I3,X2,p2,I3)) * ...
    %(kron(I3,I3,p0,I3) + kron(X1,I3,p1,I3) + kron(X2,I3,p2,I3))^2
    (kron(I3,I3,p0,I3) + kron(X1,I3,p1,I3) + kron(X2,I3,p2,I3))^2 * ...
    (kron(I3,I3,p0,I3) + kron(I3,X1,p1,I3) + kron(I3,X2,p2,I3)) * ...
    (kron(I3,I3,I3,p0) + kron(X1,I3,I3,p1) + kron(X2,I3,I3,p2)) * ...
    (kron(I3,I3,I3,p0) + kron(I3,X1,I3,p1) + kron(I3,X2,I3,p2)) * ...
    kron(I3,I3,F3,F3)*3;

    size(A);

    %ket0000 = kron(ket0, ket0, ket0, ket0);
    %A = reshape(A*ket0000,9,9)
    P = reshape(A(:,1),9,9);

end

