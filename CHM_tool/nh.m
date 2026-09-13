function n=nh(X)
% 2022-08-08

    d = size(X, 1);
    n = norm(X*X' - d*eye(d), "fro");
end
