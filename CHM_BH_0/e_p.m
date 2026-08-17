function r = e_p(U)
% 2017-10-23
% 2021-08-01 normalization to 1

    d = floor(sqrt(size(U, 1)));

    S = swap(d);
    r = SL(reshuffle(U)) + SL(reshuffle(U * S)) - 1;

%    r = r * 36/49 * 7/5; % d^2 / (d + 1)^2
end

