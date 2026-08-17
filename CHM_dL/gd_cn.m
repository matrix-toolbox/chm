function [co] = gd_cn(ci, dc)
% 2018-03-16
% 2018-03-18
% counter - to be used in gd.m (generic defect)
%
% input:
%   ci = initial counter state, eg: [0,0,0] or [0,1,0.2]
%   dc = step in (0, 1) -- should be "a proper" fraction of 1, eg. 0.1, 0.2 but not 0.3
%
% %% call:
% >> cn([0,0,0,0], 0.1)
% returns next step: [1,1,0] --> [0,0,0.5] or [-1] when overflow: [1,1,1] --> [-1]
%
% %% example of application:
%    counter=zeros(1, 4);
%    r = counter;
%    while (r ~= [-1])
%        r = cn(counter, 0.2);
%        counter = r
%    end

    cc = size(ci, 2); % number of components

    % find next non-full cell
    k = 1;
    while (abs(ci(1, k) - 1) < 1e-1 && k < cc)
        k = k + 1;
    end
    % stop, if this was the last one
    if (k == cc && abs(ci(1, k) - 1) < 1e-10)
        co = [-1];
        return
    end
    % if not - bump this (non-full) cell
    ci(1, k) = ci(1, k) + dc;
    % and reset previous cells
    for j = 1 : k - 1
        ci(1, j) = 0;
    end

    co = ci;

end
