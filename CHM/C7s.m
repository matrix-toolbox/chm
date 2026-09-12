function H=C7s
% ------------------------------------------------------------------------------
% 2022-04-23 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% 2022-04-23
%
% [1] https://arxiv.org/pdf/2204.11727.pdf
%
% Something's fishy here...
% Originally, a circulant matrix Y was obtained from sinkhorn.m, cf. ~ L7_SEED_CIRCULANT.nb script.
% Dephased form of Y was transformed into the analytic form, C', and then slightly permuted to be as close to C7D as possible.
% When C' is transformed back to a full circulant form by means of (a,b,c) --> (a, a*b, a*b*c),
% subsequent dephasing does not return the original matrix... :o
%
% Anyway, C' looks like C7D, however no transformation of equivalence is known yet.
% C' and C7D are equivalent up to a "second" core.
%
% Setting a=1 one can recover C7A or C7B...
%
% >> H = C7s;
% ------------------------------------------------------------------------------

    [zeta info] = fsolve(@p, 1.2);

    global a
    a =-1/2/zeta + 1j*( sqrt(2*sqrt(21)+6)*cos(3*pi/14)/2 - sqrt(2*sqrt(21)-6)*sin(3*pi/14) - 3* sqrt(2*sqrt(21)-6)/4);

    do
        [y info] = fsolve(@bc, exp(2j*pi*rand(1,2)));
        b = y(1);
        c = y(2);
        arg_b = mod(arg(b)+2*pi, 2*pi);
        arg_c = mod(arg(c)+2*pi, 2*pi); % get only prticular colution!
    until (abs(arg_b - 1.35) < 1e-2) && (abs(arg_c - 1.90) < 1e-2)
    do
        [y info] = fsolve(@bc, y);
        b = y(1);
        c = y(2);
        H = CHM(b, c);
    until isCHM(H)

    A=a;
    B=a*b;
    C=a*b*c;

    FCF=[ % full cyclic form
	1 A B C C B A;
	A 1 A B C C B;
	B A 1 A B C C;
	C B A 1 A B C;
	C C B A 1 A B;
	B C C B A 1 A;
	A B C C B A 1;
    ]

    %r = rand; DM=diag(exp(2j*pi*[r r r r r r r])); FCF=FCF*DM
    FCF_dephased = dephase(FCF)


end

% ------------------------------------------------------------------------------

function Y=CHM(b,c)
    global a;
    p = -a / (a+b+a*b*(1+a+b+c));
    q = (1+a+a*a/c-a/b/c+1/b/b/c-a*a/c/c-1/b/c/c+a/b/b/c/c)/(1+a);
    r = -(1+a+1/c+1/b/c+1/b/b/c)/(1+a);

    Y = [
	1 1   1     1     1     1     1;
	1 a*a a*b   a*c   a     a/c   a/b;
	1 a*b b*c   c     1/c   1/b/c 1/a/b;
	1 a*c c     1/b   a*c*r c*r   1/b/b;
	1 a   1/c   a*c*r r     a*r   1/b/c;
	1 a/c 1/b/c c*r   a*r   c*q   1/b;
	1 a/b 1/a/b 1/b/b 1/b/c 1/b   p;
    ];
end

% ------------------------------------------------------------------------------

function y=bc(x)
    global a
    % some equivalent form of the system of nonlinear equations
    y(1) = a*x(1) + a*x(2) + x(1)*x(2) + a*x(1)*x(2) + a*a*x(1)*x(2) + a*x(1)*x(1)*x(2) + a*x(1)*x(2)*x(2);
    y(2) = a + a*x(1) + x(2) + a*x(1)*x(2) + a*a*x(1)*x(1)*x(2) + a*x(1)*x(2)*x(2) + a*x(1)*x(1)*x(2)*x(2);
end

% ------------------------------------------------------------------------------

function y=p(x)
    y(1) = 1 + 12*x(1) + 18*x(1)^2 - 64*x(1)^3 - 96*x(1)^4 + 45 * x(1)^5 + 43*x(1)^6;
end

% ------------------------------------------------------------------------------

function r=isCHM(H)
    r = 0;
    if nh(H)<1e-10 && n1(H)<1e-10
        r = 1;
    end
end
