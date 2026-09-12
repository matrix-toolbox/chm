function A = A36(a=[0 0 1/3 0 0], p=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0])
% -----------------------------------------------------------------------------
% 2023-06-12
% w.bruzda[at]uj.edu.pl
% -----------------------------------------------------------------------------
% An analytic form of non-affine family of 2-unitary CHM of order 36.
% It depends on "5 + 19" complex unimodular parameters.
% In particular, when phases of the 5 first parameters are set to
%     a = [0, 0, 2*pi/3, 0, 0] % or 1/3 --- to be checked!!!
% the family reduces to an affine family depending on 19 parameters, which
% stems from the Butson matrix BH(36, 6); cf. "A36_00p00.m" and "A36_BH.m".
%
% It is based on numerical data by Suhail Ahmad Rather (2023-06-06)
% and the modification of the Sinkhorn algorithm presented in arXiv:2306.00999.
%
% Usage:
%
% >> Y = A36(a, p)
%    where a is a vector of PARTICULAR phases in R^5
%          p is a vector of ANY phases in R^19
%
% >> Y = A36([0 0 2*pi/3 0 0]/2/pi, rand(1,19));
%
% >> a(1) =  0.8320499612731626/2/pi; % some valid phases
% >> a(2) = -1.3945937333904541/2/pi;
% >> a(3) =  0.3268382320785901/2/pi;
% >> a(4) =  0.9448725457359911/2/pi;
% >> a(5) =  0.7551516958511726/2/pi;
% >> Y = A36(a, rand(1,19));
%
% An attempt to provide analytic dependencies between five additional parameters in A36 suggests they exhibit non-affine character.
% When considering a special scenario with the following pattern of fixed elements:
%
%     a = [a(1) ... a(5)] = [FIXED FIXED free free free]
%
% one always observes TWO different values for a(3) and the relation between remaining a(4) and a(5) is always linear.
% This suggests that two first parameters are additional degrees of freedom and the rest is a function of a(1) and a(2).
% -----------------------------------------------------------------------------
% special case:
% >> Y = A36([0 0 1/3 0 0], zeros(1, 19)); % is BH(36, 6)
% see also "A36_BH.m".
% -----------------------------------------------------------------------------
% 2023-06-14
% Dependencies among the phases a(j) are non-trivial and can be explored only
% numerically...
% -----------------------------------------------------------------------------

    mm = [1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1;];
    mm = [mm; shift(mm, 1)];
    minuses = []; % checkered matrix of alternating signs
    for _=1:9
        minuses = [minuses; mm];
    end

    w = exp(2j*pi/6); % SIXth roots of unity - sounds good for AME(4, SIX)
    p = exp(2j*pi*p); % change real phases into unimodular numbers
    x = exp(2j*pi*a);

    b1 = -x(1)*p(1)*x(3)/x(2);
    b2 = x(1)*p(2)*x(4)/x(5)*w^2;
    b3 = p(1)*x(3)/x(2)*w^2;
    b4 = x(1)^2*p(1)*x(3)^2/x(2)^2*w^4;
    b5 = p(2)*x(4)^2/x(5)^2*w^4;
    b6 = -x(1)*p(1)*x(3)*x(5)/x(2)/x(4)*w^2;
    b7 = p(2)*x(2)*x(4)^2/x(1)/x(3)/x(5)^2*w^2;
    b8 = p(5)/x(1)^2;
    b9 = x(2)*x(4)/x(3)/x(5);
    b10 = -x(2)/x(1)/x(3);
    b11 = b7*x(5)/x(4);
    b12 = b7*p(2)/b11;

    K1 = [p(10); p(11); p(10); p(11); p(10); p(11)];
    K2 = [p(12); p(13); p(12); p(13); p(12); p(13)];
    g0 = [p(14); p(15); p(16); p(17); p(18); p(19)];
    g1 = K1.^(1/3).*g0.*[w^2; w^2; 1; 1; w^4; w^4];
    g3 = g1./K2.^(1/3).*[w^2; w^4; 1; w^2; w^4; 1];
    G = [g1 g0 g3];


    sbB1 = zeros(18, 18);
    sbB1(7:12, 10:12) = G; % generator


    sbB1(1:6, 11:11) = sbB1(7:12, 10:10) .* [x(3)*p(7)/x(2); p(7); x(3)*p(7)/x(2); p(7); x(3)*p(7)/x(2); p(7)];
    sbB1(13:18, 12:12) = sbB1(7:12, 10:10) .* [x(4); x(5); x(4); x(5); x(4); x(5)];
    sbB1(1:6, 12:12) = sbB1(7:12, 11:11) .* [x(2); x(3); x(2); x(3); x(2); x(3)];
    sbB1(13:18, 10:10) = sbB1(7:12, 11:11) .* [p(9); x(1)^2*x(3)^2*x(5)*p(9)/x(2)^2/x(4); p(9); x(1)^2*x(3)^2*x(5)*p(9)/x(2)^2/x(4); p(9); x(1)^2*x(3)^2*x(5)*p(9)/x(2)^2/x(4)];
    sbB1(1:6, 10:10) = sbB1(7:12, 12:12) .* [x(1)^2*p(6)*x(3)/x(2); p(6); x(1)^2*p(6)*x(3)/x(2); p(6); x(1)^2*p(6)*x(3)/x(2); p(6)];
    sbB1(13:18, 11:11) = sbB1(7:12, 12:12) .* [p(8); x(2)^2*x(4)*p(8)/x(1)^2/x(3)^2/x(5); p(8); x(2)^2*x(4)*p(8)/x(1)^2/x(3)^2/x(5); p(8); x(2)^2*x(4)*p(8)/x(1)^2/x(3)^2/x(5)];

    sbB1(1:6, 7:9) = sbB1(1:6, 10:12) .* [p(3)*b8/p(5), p(4)*p(5)/b8, b8; p(3), p(4), p(5); p(3)*b8/p(5), p(4)*p(5)/b8, b8; p(3), p(4), p(5); p(3)*b8/p(5), p(4)*p(5)/b8, b8; p(3), p(4), p(5)];
    sbB1(7:12, 7:9) = sbB1(7:12, 10:12) .* [-x(3)*p(3)/x(2), -p(4)*x(2)/x(3), -p(5)*x(3)/x(2); -p(3)*b8*x(2)/x(3)/p(5), -x(3)*p(4)*p(5)/b8/x(2), -b8*x(2)/x(3); -p(3)*x(3)/x(2), -p(4)*x(2)/x(3), -p(5)*x(3)/x(2); -p(3)*b8*x(2)/x(3)/p(5), -x(3)*p(4)*p(5)/b8/x(2), -b8*x(2)/x(3); -p(3)*x(3)/x(2), -p(4)*x(2)/x(3), -p(5)*x(3)/x(2); -p(3)*b8*x(2)/x(3)/p(5), -x(3)*p(4)*p(5)/b8/x(2), -b8*x(2)/x(3)];
    sbB1(13:18, 7:9) = sbB1(13:18, 10:12) .* [x(3)*p(3)*x(5)/x(2)/x(4), x(4)*p(4)*x(2)/x(3)/x(5), x(3)*p(5)*x(5)/x(2)/x(4); x(4)*p(3)*b8*x(2)/x(3)/p(5)/x(5), x(3)*p(4)*p(5)*x(5)/x(4)/b8/x(2), x(4)*b8*x(2)/x(3)/x(5); x(3)*p(3)*x(5)/x(2)/x(4), x(4)*p(4)*x(2)/x(3)/x(5), x(3)*p(5)*x(5)/x(2)/x(4); x(4)*p(3)*b8*x(2)/x(3)/p(5)/x(5), x(3)*p(4)*p(5)*x(5)/x(4)/b8/x(2), x(4)*b8*x(2)/x(3)/x(5); x(3)*p(3)*x(5)/x(2)/x(4), x(4)*p(4)*x(2)/x(3)/x(5), x(3)*p(5)*x(5)/x(2)/x(4); x(4)*p(3)*b8*x(2)/x(3)/p(5)/x(5), x(3)*p(4)*p(5)*x(5)/x(4)/b8/x(2), x(4)*b8*x(2)/x(3)/x(5)];
    sbB1(13:18, 13:18) = sbB1(13:18, 7:12) .* [b7*w^2, p(2)*w^2, b7*w^4, p(2)*w^4, b7, p(2); p(2)*w, -b7, -p(2), b7*w^5, p(2)*w^5, -b7*w^4; b7*w^2, p(2)*w^8, b7*w^4, p(2)*w^4, b7, p(2); p(2)*w, -b7, -p(2), b7*w^5, p(2)*w^5, b7*w; b7*w^2, p(2)*w^2, b7*w^4, p(2)*w^4, b7, p(2); p(2)*w, -b7, -p(2), b7*w^5, p(2)*w^5, b7*w];
    sbB1(13:18, 1:6) = sbB1(13:18, 7:12) .* [b1*w^4, b6*w^4, b1*w^2, b6*w^2, b1, b6; b6*w^5, -b1, -b6, b1*w, b6*w, -b1*w^2; b1*w^4, b6*w^4, b1*w^2, b6*w^2, b1, b6; b6*w^5, -b1, -b6, b1*w, b6*w, -b1*w^2; b1*w^4, b6*w^4, b1*w^2, b6*w^2, b1, b6; b6*w^5, -b1, -b6, b1*w, b6*w, b1*w^5];
    sbB1(7:12, 13:18) = sbB1(7:12, 7:12) .* [b12*w^2, b5*w^2, b12*w^4, b5*w^4, b12, b5; b5*w, -b12, -b5, b12*w^5, b5*w^5, b12*w; b12*w^2, b5*w^2, b12*w^4, b5*w^4, b12, b5; b5*w, -b12, -b5, b12*w^5, b5*w^5, b12*w; b12*w^2, b5*w^2, b12*w^4, b5*w^4, b12, b5; b5*w, -b12, -b5, b12*w^5, b5*w^5, b12*w];
    sbB1(1:6, 13:18) = sbB1(1:6, 7:12) .* [b2*w^4, b11*w^2, b2, -b11*w, b2*w^2, b11; b11*w, b2*w^5, -b11, b2*w, b11*w^5, -b2; b2*w^4, b11*w^2, b2, -b11*w, b2*w^2, b11; b11*w, b2*w^5, -b11, b2*w, b11*w^5, -b2; b2*w^4, b11*w^2, b2, -b11*w, b2*w^2, b11; b11*w, b2*w^5, -b11, b2*w, b11*w^5, -b2];
    sbB1(1:6, 1:6) = sbB1(1:6, 7:12) .* [b1*w^2, p(1)*w, b1, p(1)*w^5, b1*w^4, -p(1); p(1)*w^2, b1*w, p(1), b1*w^5, p(1)*w^4, -b1; b1*w^2, p(1)*w, b1, p(1)*w^5, b1*w^4, -p(1); p(1)*w^2, b1*w, p(1), b1*w^5, p(1)*w^4, -b1; b1*w^2, p(1)*w, b1, p(1)*w^5, b1*w^4, -p(1); p(1)*w^2, b1*w, p(1), b1*w^5, p(1)*w^4, -b1];
    sbB1(7:12, 1:6) = sbB1(7:12, 7:12) .* [b3*w, b4*w, b3*w^5, b4*w^5, -b3, -b4; b4*w^2, b3, b4, b3*w^4, b4*w^4, b3*w^2; b3*w, b4*w, b3*w^5, b4*w^5, -b3, -b4; b4*w^2, b3, b4, b3*w^4, b4*w^4, b3*w^2; b3*w, b4*w, b3*w^5, b4*w^5, -b3, -b4; b4*w^2, b3, b4, b3*w^4, b4*w^4, b3*w^2];

    u1 = sbB1(1:6, 1:6);
    u2 = sbB1(1:6, 7:12);
    u3 = sbB1(1:6, 13:18);
    u4 = sbB1(7:12, 1:6);
    u5 = sbB1(7:12, 7:12);
    u6 = sbB1(7:12, 13:18);
    u7 = sbB1(13:18, 1:6);
    u8 = sbB1(13:18, 7:12);
    u9 = sbB1(13:18, 13:18);

    dd1 = [x(1), x(4)/x(5), x(1), x(4)/x(5), x(1), x(4)/x(5); -x(4)/x(5), -x(1), -x(4)/x(5), -x(1), -x(4)/x(5), -x(1); x(1), x(4)/x(5), x(1), x(4)/x(5), x(1), x(4)/x(5); -x(4)/x(5), -x(1), -x(4)/x(5), -x(1), -x(4)/x(5), -x(1); x(1), x(4)/x(5), x(1), x(4)/x(5), x(1), x(4)/x(5); -x(4)/x(5), -x(1), -x(4)/x(5), -x(1), -x(4)/x(5), -x(1)];
    dd2 = [b9, -b10, b9, -b10, b9, -b10; b10, -b9, b10, -b9, b10, -b9; b9, -b10, b9, -b10, b9, -b10; b10, -b9, b10, -b9, b10, -b9; b9, -b10, b9, -b10, b9, -b10; b10, -b9, b10, -b9, b10, -b9];
    dd3 = [x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2), x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2), x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2); -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5), -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5), -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5); x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2), x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2), x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2); -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5), -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5), -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5); x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2), x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2), x(2)*x(4)/x(1)/x(3)/x(5), x(3)/x(2); -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5), -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5), -x(3)/x(2), -x(2)*x(4)/x(1)/x(3)/x(5)];

    sbA1 = [dd1.*u1, dd2.*u2, dd3.*u3;     dd2.*u4, dd3.*u5, conj(dd1).*u6;     dd3.*u7, conj(dd1).*u8, conj(dd2).*u9];
    sbA2 = sbA1.*minuses;
    sbB2 = sbB1.*minuses;

    A = [sbA1, sbB1; sbB2, sbA2];

end


    % clf
    % colormap(prism(6))
    % imagesc(LF(Y))
    % axis square
    % axis off
    %
    % additional script is required to put numbers on colored cells



%    hint for a future:
%
%    domain = 2*pi;  % for a certain reason it is better to work in [0..2pi) interval, not [0..1)
%                    % in that case "mu" can be as large as 0.01 and convergence is much faster
%                    % however, one must remember about rescaling phases "a" when being passed to A36!
%                    % ALSO fixed entries (in f5) must (should) belong to [0..2pi).



