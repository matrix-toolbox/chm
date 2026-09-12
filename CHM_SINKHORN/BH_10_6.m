function Y = BH_10_6()
% ------------------------------------------------------------------------------
% 2023-02-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Hermitian BH(10, 6) with generic d = 8 and #L = 6.
% Found as a solution of "HH_10_12_1002.m" and defined as a separate pattern.
% ------------------------------------------------------------------------------

    s3 = sqrt(3) / 2;

    a =  0.5 - s3*i;
    b = -1.0;
    c =  0.5 - s3*i;
    d =  0.5 + s3*i;
    e =  0.5 + s3*i;
    f =  0.5 + s3*i;
    g = -1.0;
    h = -0.5 + s3*i;
    j =  0.5 + s3*i;
    k = -0.5 + s3*i;
    m = -0.5 - s3*i;
    n =  1.0;
    o =  0.5 - s3*i;
    p = -1.0;
    q =  0.5 + s3*i;
    r =  0.5 - s3*i;
    s =  0.5 - s3*i;
    t =  0.5 + s3*i;
    u =  0.5 - s3*i;
    v =  0.5 + s3*i;
    w =  0.5 + s3*i;
    x =  0.5 - s3*i;
    y =  0.5 - s3*i;
    z = -1.0;

    Y = [
        1  1    1    1    1    1    1      1      1      1   ;
        1  1   -1   -1   -1    a    b      c      d      e   ;
        1 -1    1   -1   -1    f    o      p      q      r   ;
        1 -1   -1    1   -1    g    s      t      u      v   ;
        1 -1   -1   -1    1    j    w      x      y      z   ;
        1  a'   f'   g'   j'  -1    h      k      m      n   ;
        1  b'   o'   s'   w'   h'  -1      k/h    m/h    n/h ;
        1  c'   p'   t'   x'   k'   h/k   -1      m/k    n/k ;
        1  d'   q'   u'   y'   m'   h/m    k/m   -1      n/m ;
        1  e'   r'   v'   z'   n'   h/n    k/n    m/n   -1   ;
    ];

end
