function n=n1(X)
% 2022-08-08

    d = size(X, 1);
    n = norm(abs(X) - ones(d), "fro");
end
