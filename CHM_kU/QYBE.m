function r = QYBE(H, display_mode)
% ------------------------------------------------------------------------------
% 2024-07-12
% 2024-07-13
%
% based on:
% https://arxiv.org/pdf/2202.12306 --> construction of dual unitary U from CHM
% https://arxiv.org/pdf/2406.03781 --> QYBE for symmetric CHM
%
% more on QYBE:
% https://arxiv.org/pdf/1704.03558
% ------------------------------------------------------------------------------
% OBSERVATIONS for some SYMMETRIC CHM
%
%	FULFILLING QYBE:
%	1) F4([0.75])
%	2) F9([8, 7, 7, 5]/9)
%	the above cases are interesting, because they fulfill QYBE "twice":
%	a) one needs not create the large U gate, becase their size is a square
%	   and this works only for these particular phases!
%	b) via U it works as well for any phase for F4 and this particular one for F9... (?)
%	
%	TODO: check other Fq(p)
%	      they all work for Fq(0), but what about p != 0?
%	
%	3) BH_6_4_6 <-- it is equivalent to F6
%	4) F6([0, 0])
%	5) F2 (x) F4(p) for any p in [0, 2pi)
%	6) F4(p1) (x) F4(p2) for any p1, p2 in [0, 2pi)
%	7) F2 (x) BH_6_4_6
%	8) F2 (x) BH_6_4_6 (x) F2 and probably for other such tensor products where each component is in QYBE class
%	9) BH(12, 6, 1)
%	10) H16A -- it is also self-G-dual
%	11) P8([0 0 0])
%
%	FAIL:
%	1) BH_6_4_4
%	2) S6 (spectral CHM)
%	3) C6 (C-6-R)
%	:
%	many others... IN PROGRESS...
% ------------------------------------------------------------------------------

    q = size(H, 1); % local dimension

    U = zeros(q^2, q^2); % unitary cell U in U(q^2) for the Hadamard circuit
    for a = 1:q
    for b = 1:q
    for c = 1:q
    for d = 1:q
        AB = a+(b-1)*q; % multi-index for U
        CD = c+(d-1)*q; % multi-index for U
        U(AB, CD) = H(a, b) * H(b, d)' * H(d, c) * H(c, a)' / q;
    end, end, end, end


    EPSILON = 1e-14;
    assert(isU(U), 1, EPSILON); % assert that U is a unitary of order q^2 up to EPSILON


    LHS = kron(U, eye(q)) * kron(eye(q), U) * kron(U, eye(q)); % check QYBE
    RHS = kron(eye(q), U) * kron(U, eye(q)) * kron(eye(q), U); % check QYBE


    r = norm(LHS - RHS, "fro");

    if strcmp(display_mode, "VERBOSE")
        message = "DOES NOT FULFILL";
        if r < 1e-8
            message = "FULFILLS";
        end
        printf("numerically it seems that this matrix %s QYBE\n", message);
    end

end

