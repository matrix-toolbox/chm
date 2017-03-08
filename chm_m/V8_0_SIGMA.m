% 20160606
% W. Bruzda, name[at]alumni.uj.edu.pl | name = w.bruzda
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/

% Eight numerical solutions of non-Butson type.
% Each row of the SIGMA_ARRAY has [a, b, c] structure.

% Courtesy of Veit Elser. Unpublished notes, 2011.

% >> SIGMA = 1 % 2, 3 or 4
% >> H = V8_0_SIGMA(SIGMA)
% >> abs(H .* H'), norm(H * H' - 8 * eye(8), 'fro')

function H = V8_0_SIGMA(SIGMA)

    SIGMA_ARRAY = [
        -0.6508910153042045364155285204797 + i * 0.75917118372358007721802306013268, ...
        -0.7799568409991842327822652885197 - i * 0.62583330542451417655381771182110, ...
        +0.6183387264655547438069198981350 + i * 0.78591171218716156084627407300094;
        -0.6508910153042045364155285204797 - i * 0.75917118372358007721802306013268, ...
        -0.7799568409991842327822652885197 + i * 0.62583330542451417655381771182110, ...
        +0.6183387264655547438069198981350 - i * 0.78591171218716156084627407300094; % cc

        +0.9741272443714635143662193012491 + i * 0.22600024728583596147147890367920, ...
        -0.6183387264655547438069198981350 - i * 0.78591171218716156084627407300094, ...
        -0.1941704033722792815839540127294 + i * 0.98096781519795355342764679247832;
        +0.9741272443714635143662193012491 - i * 0.22600024728583596147147890367920, ...
        -0.6183387264655547438069198981350 + i * 0.78591171218716156084627407300094, ...
        -0.1941704033722792815839540127294 - i * 0.98096781519795355342764679247832; % cc

        -0.9741272443714635143662193012491 - i * 0.22600024728583596147147890367920, ...
        +0.0325522888386497926086086223447 - i * 0.99947003381360319608224059999136, ...
        +0.7799568409991842327822652885197 - i * 0.62583330542451417655381771182110;
        -0.9741272443714635143662193012491 + i * 0.22600024728583596147147890367920, ...
        +0.0325522888386497926086086223447 + i * 0.99947003381360319608224059999136, ...
        +0.7799568409991842327822652885197 + i * 0.62583330542451417655381771182110; % cc

        +0.6508910153042045364155285204797 + i * 0.75917118372358007721802306013268, ...
        +0.1941704033722792815839540127294 - i * 0.98096781519795355342764679247832, ...
        -0.0325522888386497926086086223447 - i * 0.99947003381360319608224059999136;
        +0.6508910153042045364155285204797 - i * 0.75917118372358007721802306013268, ...
        +0.1941704033722792815839540127294 + i * 0.98096781519795355342764679247832, ...
        -0.0325522888386497926086086223447 + i * 0.99947003381360319608224059999136; % cc
    ];

    row = zeros(1, 3);
    try
        row = SIGMA_ARRAY(SIGMA, :);
    catch
        warning('No valid SIGMA value provided! Continue with default SIGMA = 1.');
        row = SIGMA_ARRAY(1, :);
    end

    a = row(1);
    b = row(2);
    c = row(3);

    V = [
        -1, -1,  b,  b,  c,  c,  a,  a;
        -1,  b, -1,  c,  b,  a,  c, -a;
         b, -1,  c, -1,  a,  b, -a,  c;
         b,  c, -1,  a, -1, -a,  b, -c;
         c,  b,  a, -1, -a, -1, -c,  b;
         c,  a,  b, -a, -1, -c, -1, -b;
         a,  c, -a,  b, -c, -1, -b, -1;
         a, -a,  c, -c,  b, -b, -1,  1;
    ];

    LD = diag([-1, -1, 1/b, 1/b, 1/c, 1/c, 1/a, 1/a]);
    RD = diag([1, 1, -1/b, -1/b, -1/c, -1/c, -1/a, -1/a]);

    H = LD * V * RD; % dephased form

end

