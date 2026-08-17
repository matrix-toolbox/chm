function Y = Y_10_3_278(p)
% -----------------------------------------------------------------------------
% 2022-08-06 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% -----------------------------------------------------------------------------
% A 3-parametric non-affine family of CHM of order N = 10 found by the Sinkhorn algorithm.
% It has generic defect = 3 and #L = 278.
% Phases p1, p2, p3 do not run over full circle!
%
% [1] https://arxiv.org/abs/2204.11727
% -----------------------------------------------------------------------------

    % free parameters:
    a = exp(2j*pi * p(1));
    b = exp(2j*pi * p(2));
    u = (-a + 2* a^2 - b + 4*a*b - a^2*b + 2*b^2 - ...
        a*b^2 - sqrt(-4*(-a - b + a*b)*(a*b - a^2*b - a*b^2) + (a - 2*a^2 +  ...
        b - 4* a* b + a^2* b - 2* b^2 + a* b^2)^2))/(2*(-a - b + a *b));


    % free parameter: (full angle)
    c = exp(2j * pi * p(3));


    % functions:
    p = - 2 - a/b - b/a;
    q = (p/2 - sqrt(p^2 - 4)/2 + 1) / (1 - 1/u - 1/a - 1/b); % here hides another branch of the root!
    r = q - 1 - q/u - q/a - q/b;
    s = (q + q * r - r) / (1 + r - q);
    t = - r * u * (a + b) / (r * u + r + s + r * s) / b;


    Y = [
        1   1     1     1       1       1         1         1       1         1;
        1   1     u    -s*q/r  -s*q/r   s*q/r     s*q/r    -u      -1        -1;
        1   1    -1     s       s/r    -s*q/b/r  -s*q/a/r   u      -1/c       1/c;
        1   1    -1     s/r     s      -s*q/a/r  -s*q/b/r   u       1/c      -1/c;
        1   r     q    -q      -q       q*c      -q*c       r*u/s   r*a/s     r*b/s;
        1   1/r   q/r  -q/r    -q/r    -q*c/r     q*c/r     u/s     b/s       a/s;
        1   a/b   a     a*c    -a*c    -a        -a         t      -r*a/s    -a/s;
        1   b/a   b    -b*c     b*c    -b        -b         t*b/a  -b/s      -r*b/s;
        1  -1    -c    -a*c    -b*c    -q*c      -q*c/r    -u*c     c         c;
        1  -1     c     b*c     a*c     q*c/r     q*c       u*c    -c        -c;
    ];


end
