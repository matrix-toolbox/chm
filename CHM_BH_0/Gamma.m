function Y = Gamma(X, ss, ssd)
% 2006
%
% Partial Transpose (c) Toby Cubitt.
%
% This code is entirely based on Toby's work.
% usage:
%
% >> d = 4;
% >> X = randn(d*d);
% >> Y = Gamma(X, 2, [d d]);

    switch length(ssd)
    case 0
        d1 = 2;
        d2 = 1;
        d3 = 2;
    case 1
        d1 = ssd(1);
        d2 = 1;
        d3 = d1;
    case 2
        d1 = ssd(1);
        d2 = 1;
        d3 = ssd(2);
    case 3
        d1 = ssd(1);
        d2 = ssd(2);
        d3 = ssd(3);
    end

    if (d2 == 1 & ss == 2) ss = 3; end

    switch ss
    case 1
        Y = X;
        i = (1 : d2 * d3);
        for j = 0 : d1 - 1
            for k = 0 : d1 - 1
                if j == k continue; end
                Y(i + j * d2 * d3, i + k * d2 * d3) = ...
                X(i + k * d2 * d3, i + j * d2 * d3);
            end
        end
    case 2
        Y = X;
        i = kron(ones(1, d1), [1: d3]) + ...
        kron(dim2 * d3 * [0 : d1 - 1], ones(1, d3));
        for j = 0 : d2 - 1
            for k = 0 : d2 - 1
                if j == k continue; end
                Y(i + j * d3, i + k * d3) = ...
                X(i + k * d3, i + j * d3);
            end
        end
    case 3
        Y = X;
        i = (1: d3);
        for j = 0 : d1 * d2 - 1
            for k = 0 : d1 * d2 - 1
                Y(i + j * d2 * d3, i + k * d2 * d3) = ...
                X(i + j * d2 * d3, i + k * d2 * d3)';
            end
        end
    end
end

