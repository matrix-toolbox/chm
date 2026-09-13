function verdict = cSH(X)
% 2023-02-03
%
% Check if X is symmetric or Hermitian.
    verdict = "";

    if norm(X - transpose(X), "fro") < 1e-7
        verdict = "matrix is symmetric: X = X^T";
    end

    if norm(X - X', "fro") < 1e-7
        verdict = "matrix is Hermitian: X = X^*";
    end

end
