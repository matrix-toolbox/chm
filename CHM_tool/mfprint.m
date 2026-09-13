function PHI =  mfprint(H)
% 2017-12-19
% matrix fingerprint
% cf. F. Szollosi, Exotic CHM and their equivalence; https://arxiv.org/pdf/1001.3062.pdf
%
% call:
% >> PHI = mfprint(H);
% >> PHI{2}
% >> PHI{3} % etc...
%
%
% 2022-04-22
% Improved comparison of numerical values.
% Now the result is more credible.

  numericalTolerance = 1e-7;

  N = size(H, 1);
  PHI = {};

  for d = 2 : floor(N / 2)
    vi = [];
    m = nchoosek([1:N], d);                 % prepare all possible minor coordinates of size d
    s1 = size(m, 1);
    s2 = size(m, 2);
    for rs = 1 : s1;
      for cs = 1 : s1;
        [a b] = meshgrid(m(rs,:), m(cs,:)); % Matlab-Cartesian product
        rc = [b(:) a(:)];                   % rc = [row column] of matrix mu = minor
        mu = zeros(s2, s2);                 % cardinality of rc should be s2^2
        mu = reshape(H(sub2ind(size(H), rc(:,1), rc(:,2))), s2, s2);
                                            % mu emerges out of H by means of rc
        vi = [vi, abs(det(mu))];            % collection of all abs(minor)
      end
    end

    vi = sort(vi);
%    ui = unique(vi)
    ui = getUnique(transpose(vi),numericalTolerance);
    ui = sort(ui);
    mi = zeros(size(ui));
    for i = 1 : length(ui)
%      mi(i) = sum(vi == ui(i));             % get multiplicities m_i(d) of each v_i(d) 
      mi(i) = sum(abs(vi - ui(i))<numericalTolerance);             % get multiplicities m_i(d) of each v_i(d) 
    end

    disp(sprintf(' v_i(%d)  m_i(%d)', d, d));
    for j=1:size(mi,2)
      disp(sprintf(' %1.14g \t %d', ui(j), mi(j)));
    end
    disp(sprintf('\n'));
    % disp(transpose([ui; mi]));

    PHI{d} = transpose([ui; mi]);           % return
  end
end

