function H = SH_9_0_76_pattern()
% ------------------------------------------------------------------------------
% 2025-02-21
% general form of a matrix ss_9_0_L76*
% ------------------------------------------------------------------------------

    addpath ../matrix_tool
    addpath ../CHM

    npa = 4; % number of parameters in a matrix

    while 1
        pa = 2*pi*rand(1,npa);
        mu = 1;
        Z_OPTIMAL = +Inf;
        while (Z_OPTIMAL > 5e-14)
            RESTORE_pa = pa;
            pa = mod(pa + mu * randn(1,npa), 2*pi);

            b = exp(1j*pa(1,1));
            c = exp(1j*pa(1,2));
            d = exp(1j*pa(1,3));
            e = exp(1j*pa(1,4));

% choose one matrix
% analytical values for parameters are currently out of reach!

%            x = -(2+2*d+e+e/c+2*d/b);
%            H = [
%                1       1       1       1        1       1         1           1        1       ;
%                1       1       1       b        c       1/b       1/c         1/b      b       ;
%                1       1       x       d        e       d/b       e/c         d/b      d       ;
%                1       b       d      -b*e      b*e    -e^2/c     e/c        -e/c     -e^2*b/c ;
%                1       c       e       b*e      e^2     e         e^2/c       e/b      e       ;
%                1       1/b     d/b    -e^2/c    e      -e/c/b     e/c/b      -e^2/b/c -e       ;
%                1       1/c     e/c     e/c      e^2/c   e/c/b     e^2/c^2     e/c      e*b/c   ;
%                1       1/b     d/b    -e/c      e/b    -e^2/b/c   e/c        -e/b     -e^2/c   ;
%                1       b       d      -e^2*b/c  e      -e         e*b/c      -e^2/c   -e*b/c   ;
%            ];

            x = -(2  + c + 2*d/e - d/e/c +  2*c/e);
            H=[
                1       1       1             1             1        1          1        1        1       ;
                1       c^2    -c^2          -1             c        d          e        c^2/d    c^2/e   ;
                1      -c^2     d^2/e^2      -d^2/c^2/e^2   d/e      d^2/e^2    1        d/e      d/e     ;
                1      -1      -d^2/c^2/e^2   d^2/c^2/e^2  -d/e/c   -d^2/e/c^2 -d/c^2   -d/e^2   -1/e     ;
                1       c       d/e          -d/e/c         x        d/e        1        c/e      c/e     ;
                1       d       d^2/e^2      -d^2/e/c^2     d/e      d/e        d/e      d/e^2   -c^2/e   ;
                1       e       1            -d/c^2         1        d/e        e/d     -c^2/d    1/e     ;
                1       c^2/d   d/e          -d/e^2         c/e      d/e^2     -c^2/d    c^2/e/d  c^2/e^2 ;
                1       c^2/e   d/e          -1/e           c/e     -c^2/e      1/e      c^2/e^2  c^2/e/d ;
           ];


            Z = nh(H);
            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z_OPTIMAL * 0.02;
            else
                pa = RESTORE_pa;
            end
        end % while 1...

        mod(pa+2*pi,2*pi)/2/pi
        return % <-------------- uncomment to get the matrix

    end
end


