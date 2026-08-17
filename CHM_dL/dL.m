function dL(fh, np, nMax)
% 2023-02-03
%
% Randomly collect values of defect and #L for a given CHM.
%  {  fh = function handler
%  {  np = #parameters that fh depend on
%  {  nMax = maximum number of phase numerator, see below
%
% Usage:
% >> dL(@F12A, 9, 24); % watch @ before argument!
%
% This is not a comprehensive search.
% It is possible that many values of d and #L were not encountered.
% When checking a matrix (usually Fouriers and other affine families),
% its parametrization is assumed to depend on phases in [0, 1).
% Phases take the form n/d with n<d.
% For special matrices, one must adjust the scope of parameters.

    addpath ../CHM
    addpath ../matrix_tool

%  if ~exist('defect')
%      error('defect is undefined! please locate ''defect'' scripts to continue!');
%  end


    KNOWN = [ % d, #L
    ];

    d = zeros(1, np);
    n = zeros(1, np);

    while 1;
        p = [];
        for j=1:np
            d(j) = randi(nMax);
            n(j) = randi(d(j))-1;
            p = [p, n(j)/d(j)];
        end;
        X = fh(p);
        L = cL(X, 1e-8);
        defect = ud(X, "S", 1e-8);

        f = 0;
        for j=1:size(KNOWN, 1) % quite inefficient loop, but it is enough...
            k = KNOWN(j, :);
            if (k(1) == defect) && (k(2) == L)
                f = 1;
                break;
            end
        end

%%%        if L<20 && !f
        if !f
            for j=1:np
                g = gcd(n(j), d(j));
                n(j) /= g;
                d(j) /= g;
                if n(j)==0
                    printf("0\t");
                else
                    printf("%d/%d\t", n(j), d(j));
                end
            end
            printf("%d\t%d\tq=%d\n", defect, L, isBH(X, 10000, ""));
            KNOWN = [KNOWN; defect, L];
        end
    end

end


% alt.
%        gd = ud(D10X_7(1, rand(1, 7)), "S", 1e-8);
%        if !any(gdD10X_7(:) == gd)
%            gdD10X_7 = [gdD10X_7 gd]
%        end

