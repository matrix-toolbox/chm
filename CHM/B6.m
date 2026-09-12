function H = B6(p)
% ------------------------------------------------------------------------------
% 2006-12-05 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% First non-affine family of CHM of order N = 6 found by K. Beauchamp and R. Nicoara;
% [1] https://arxiv.org/abs/math/0609076
%
% Watch the scope of p: a < p < 1 - a for a = acos((sqrt(3) - 1) / 2) / 2/pi;
%
% >> H = B6(1/3);
% ------------------------------------------------------------------------------

    a = acos((sqrt(3) - 1) / 2) / 2/pi;
    if p < a || p > 1-a
        error("Invalid parameter!");
    end

    y = exp(2j*pi*p);

    s = 2*randi(2) - 3; % choose either branch...
    x = (1+2*y+y^2 - s*sqrt(2)*sqrt(1+2*y+2*y^3+y^4)) / (1+2*y-y^2);
    z = (1+2*y-y^2)/y/(y^2+2*y-1);
    t = x*y*z;

    H = [
        1   1   1   1   1   1  ;
        1  -1  -x' -y   y   x' ;
        1  -x   1   y   z' -t' ;
        1  -y'  y' -1  -t'  t' ;
        1   y'  z  -t   1  -x' ;
        1   x  -t   t  -x  -1  ;
    ];

%    OLD VERSION - ?
%
%    q1 = acos((sqrt(3) - 1) / 2) / (2 * pi); % = 0.19035916268767
%    q2 = 1 - acos((sqrt(3) - 1) / 2) / (2 * pi); % = 0.80964083731233
%    t = exp(2j * pi * (q1 + p * (q2 - q1)));
%
%    Bn = (1 + 2*t + t*t - sqrt(2)*sqrt(1 + 2*t + 2*t^3 + t^4));
%    Bx = Bn / (2*t + 1 - t*t);
%    Bz = Bn / (2*t - 1 + t*t);
%    By = (2*t + 1 - t*t) / (t * (2*t - 1 + t*t));
%
%    H = [
%        1  1   1    1   1    1  ;
%        1 -1  -Bx' -t   t    Bx';
%        1 -Bx  1    t   By' -Bz';
%        1 -t'  t'  -1  -Bz'  Bz';
%        1  t'  By  -Bz  1   -Bx';
%        1  Bx -Bz   Bz -Bx  -1  ;
%    ];

end
