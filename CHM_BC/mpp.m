function n0 = mpp(pc = 0, pd = 0, a1)
% ------------------------------------------------------------------------------
% 2025-06-16
%
% calculate roots of a monic palindromic polynomial (mpp)
% and compare them with entries of a given (isolated) CHM
% (actually phases are compared)
%
% >> mpp(pc = 0, pd = 0, a1)
%
% pd should be > 1 to get non-trivial results
% pc, pd are OPTIONAL while a1 is MANDATORY (a1 = phases of unimodular entries)
%
% >> load ../CHM_symmetric/ss_11_0_139_....data
% >> N = size(A, 1); A = reshape(A(2:N, 2:N), (N-1)^2, 1); % vertical vector
% >> n0 = 0; while n0==0, n0 = mpp(0, 0, A); end
%
% ------------------------------------------------------------------------------


    if pc == 0 % no coefficients ==> creatge random palindromic polynomial (pp)
        if !pd, pd = 1 + randi(15); else pd++; end
%       printf("random mpp of (random) degree: %d\n", pd);
        pc = rpc(pd + 1, -4096, 4096);
    else % create pp from pc (pd is ignored)
        if !check_p(pc), error("it is not a palindromic sequence!\n"); end
    end

    rp = roots(pc);
    r1 = @(r, Epsilon) r(abs(abs(r) - 1.0) < (nargin < 2) * 1e-8 + (nargin == 2) * Epsilon); % select only unimodular roots
    p1 = mod(arg(r1(rp, 1e-10)) + 2 * pi, 2 * pi) / 2 / pi; % get their phases and ...
    n0 = any(any(abs(p1(:) - a1(:)') < 1e-12)); % ... check if there is a non-empty intersection of p1 and a1 up to 1e-12
    if n0, pc, beep2(); end % if so, show pc and RETURN 

end

function c = rpc(d, m0, m1)
    % generates a random palindromic sequence of length 1 < d <= pd in [m0, m1]
    % with first (and last) element equal to 1 (monic pp)
    L = ceil(d / 2);
    r = randi([m0, m1], 1, L);
    r(1) = 1; % monic
    c = [r, fliplr(r(1 : floor(d / 2)))];
end

function r = check_p(ppc)
    % check if the sequence ppc is palindromic
    r = isequal(ppc, flip(ppc));
end





% # apt install octave-dev
% # octave
% octave:> pkg install -forge struct
% octave:> pkg install -forge statistics
	% throws warnings like:
	% warning: function /home/x/.local/share/octave/api-v58/packages/statistics-1.7.4/shadow9/var.m shadows a core library function
% octave:> pkg install -forge optim
% octave:> pkg load optim



