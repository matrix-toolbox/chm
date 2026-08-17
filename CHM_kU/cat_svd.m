function cat_svd(d, a, b, c)
% ------------------------------------------------------------------------------
% 2022-12-27 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Quantized cat map - Singular Value Decomposition.
% ------------------------------------------------------------------------------

    N = d^2;

    F = cat(d,a,b,c);

    U = F*flipud(eye(N))/d; % == fliplr(F)/d; reverse columns of F
    D = d*eye(N);
    V = flipud(eye(N));

    F;
    U*D*V';
    SVD_CORRECT = norm(F-U*D*V',"fro")

    norm_h = nh(F)
    norm_1 = n1(F)
    SL3(F)

end

