function H = TU_16_19_26(p)
% ------------------------------------------------------------------------------
% 2022-11-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Two-parametric family which yields 2-unitary CHM for particular values
% of p1 and p2.
%
% Originally obtained value producing 2-unitary CHM:
%    b =    0.04836146778001216 +     0.99882989964966758*1j;   % arg(b) = 0.2423000298673765
%    d =   -0.02737385355500471 +     0.99962526585793490*1j;   % arg(d) = 0.2543572283857645
%
% TO BE CONTINUED...
% ------------------------------------------------------------------------------

    b = exp(2j*pi*p(1));
    d = exp(2j*pi*p(2));

    H = [
        1       1       1       1       1       1       1       1         1       1       1       1        1       1      1       1      ;
        1       1       b       b^2    -1      -b      -b^2    -b^2       b^2     b      -b      -b        b       b     -b      -1      ;
        1       1/b     b       b      -1/b    -1/b    -b      -1        -1      -1       1/b     1       -1      -b      1       1      ;
        1      -1      -b^2     b^2     1      -1      -b^2     b^2      -b       b       b      -b       -b       b      b      -b      ;
        1       d      -b*d    -b       d       1      -b      -b*d       b*d    -d      -1       b       -d       b*d    b      -1      ;
        1      -1      -b      -b*d    -1      -d      -b*d     b*d       b      -1       1       b*d      1       b*d   -b*d     d      ;
        1       1/d    -b       b       1/d    -1/d     b       b/d      -b      -1/d    -1       b        1/d    -b/d   -b      -1/d    ;
        1      -1/d     b/d    -b       1/d    -1       b      -b/d       b       1       1/d     b/d     -1      -b     -b/d    -1/d    ;
        1       1      -1      -1       1       1      -1      -1        -1       1       1      -1        1      -1     -1       1      ;
        1       1       b      -b^2    -1       b       b^2     b^2      -b^2    -b       b      -b       -b       b     -b      -1      ;
        1      -1/b     b       b       1/b     1/b    -b       1         1      -1      -1/b    -1       -1      -b     -1       1      ;
        1      -1       b^2    -b^2     1      -1       b^2    -b^2      -b      -b      -b      -b        b       b      b       b      ;
        1      -d       b*d    -b      -d       1      -b       b*d      -b*d     d      -1       b        d      -b*d    b      -1      ;
        1      -1      -b       b*d    -1       d       b*d    -b*d       b      -1       1      -b*d      1      -b*d    b*d    -d      ;
        1      -1/d    -b       b      -1/d     1/d     b      -b/d      -b       1/d    -1       b       -1/d     b/d   -b       1/d    ;
        1       1/d    -b/d    -b      -1/d    -1       b       b/d       b       1      -1/d    -b/d     -1      -b      b/d     1/d    ;
    ];

end

