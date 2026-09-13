function A = QYBE_RW(Mh, pn)
% ------------------------------------------------------------------------------
% 2024-07-16
%
% Cf.
% [1] papers by J. Hietarinta
% [2] arXiv: 2407.10731
% [3] papers by Agata and Alicja Smoktunowicz (informed by W. Tadej)
%
%
% It seems that in the CHM class, QYBE holds ONLY for matrices equivalent
% to the Fourier one. This explains the case that BH_6_4_4 does not work while
% the other Butson-like, BH_6_4_6 ~ F_6, does.
% This fact should be easily proved and it should be written somewhere
% in the paper by PWC and AL.
%
%
% ------------------------------------------------------------------------------
% input:
%     Mh = matrix handler for M
%     pn = numer of affine parameters in M
% ------------------------------------------------------------------------------
% >> QYBE_RW(@F4, 1) --> any phase A in [0, 1) is "optimal"
% ------------------------------------------------------------------------------
% >> QYBE_RW(@F6, 2) --> optimal phases A = [0 0]
% >> @M6(A) for A = 1
% >> @D6(A) for no A...
% >> @K6_2(A) for no A...
% >> @K6_3(A) for no A...
% >> @B6(A) for no A...
% ------------------------------------------------------------------------------
% >> @P7(1) does not work at all
% ------------------------------------------------------------------------------
% >> QYBE_RW(@F8, 5) --> optimal phases:
%    1)
%    0.09234009784141135   1/4   0.09234009784141135+3/4   1/4   0.09234009784141135+3/4
%    0.09234009784141135 <-- should be recovered analytically
%
%    2)
%    0.8781803725441183    1/2   A(1)-1/2                  1/2   A(1)-1/2
%
%    it is very likely that these paremeters make the Fourier symmetric
% ------------------------------------------------------------------------------
% >> QYBE_RW(@F10, 4) --> optimal phases:
%    1)
%    trivial ones: 0 1 1 0
%                  0 0 1 1
%                  1 1 1 1
% ------------------------------------------------------------------------------
% >> QYBE_RW(@F12B, 9)
%    ...
%
% ------------------------------------------------------------------------------

 

    Z_OPTIMAL = +Inf;
    mu = 1;
    A = rand(1, pn); % random initial phases for M = M(A)
    H = Mh;

    while Z_OPTIMAL > 2e-13

        RESTORE_A = A;
        A = mod(A + mu*randn(1, pn), 1);
        Z = QYBE(H(A), "");

        if Z < Z_OPTIMAL
            Z_OPTIMAL = Z
            mu = Z * 0.001;
        else
            A = RESTORE_A;
        end

    end


end
