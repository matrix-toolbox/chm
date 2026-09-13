function H = dita_construction(A, varargin)
% -------------------------------------------------------------------------------
% 2026-08-13 m$c
% -------------------------------------------------------------------------------
  K = size(A, 1);

  if ndims(A) != 2 || size(A, 2) != K
    error("A must be a square matrix");
  endif

  if numel(varargin) != K
    error("Exactly size(A,1) matrices B_1,...,B_K are required");
  endif

  M = size(varargin{1}, 1);

  for j = 1:K
    if ndims(varargin{j}) != 2 || rows(varargin{j}) != M || columns(varargin{j}) != M
      error("All B_j matrices must be square and have the same size");
    endif
  endfor

  blocks = cell(K, K);

  for i = 1:K
    for j = 1:K
      blocks{i, j} = A(i, j) * varargin{j};
    endfor
  endfor

  H = cell2mat(blocks);
endfunction
