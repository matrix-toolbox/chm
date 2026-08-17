function H=sinkhornY11(kMax)
% 2022-11-08
%
% sinkhornY11(8192)
%
%
% the bigger kMax the more findings...
%    octave:2> sinkhornY11(2^16)
%    xxxxxxxx.xxxxx.xxx.xxxxxxxx...xxx.xx....xx.xxxxxxx.xxxxxx..xxxxxxx.xxxxxx.xx.x.xxxx.xxx.xx.xx..xxx..xx..x.xxxxxxxx.xx.xxx.x.xxx..xxx..xxxxxxxxxxxxxxxxxx.xx.x.xxx..xxx...x.xxx.xx.xxxxxxx.xx...xxxxxxxxxxxxxxx..xxxx.xxxxx.xxxxxx...x.xx.xxxxxx.xxx.xxxx....x.xxxxx..xxxxxx..xx..xx.xxxxxxxxxx.x.xxx.xxx.x.xxxxxxx.x.x.xxxxx.xx..xxxx.x.x...x.x.x..xxxx.xxx.x..xxx..x.x.xxxxxxxx..x.xxxxx.xxx.xxx.xxxx.xx.x.x.xxx.x...xx.x.xx.xxxxxxx.xxxxxx.xxxx.xxxx.xx^C
%    octave:3> sinkhornY11(2^17)
%    x...xxxxxxxxx.xxxxx...xxxxxxxxx.x..xx....x.x.xxx.xxx.x.xxx.x...x..x......x.xx.x..xx....xx..x..xxx.xxxx.xx.xxx.x.xxxxxx....xx..x.x......xxxx..x...x..xx....x..x.xx..xx..xxx.xxxxx...xx..x...xx...xxx.xx.x..xxx.xxxx.xx.x...x....xxxxxxxx..xxx..x.x.xxx.x.xxx..x.xxx.xx....xxx..x.xx.x...x.xxxx.xx...xxxx.xxxx.....xx.xx..x....x..x.x.x.xxx...xxxx.xx.xx.xx.x..xxxxx..x.xxx.x....x.xx.xxxxxx..x.xxxxxxxx.x.xx..xxxxx....xxx...x...xxx...x..xx.xxx.xxx..x.xxxx.xxx..xxx....x.xx.xx..xxxx^C^C
%    octave:6> sinkhornY11(2^18)
%    x.......
%    octave:7> sinkhornY11(2^19)
%    x..x..x..xx..x..x.x..x.x.xx..x.x..xx
%    saved as Y11_0_L323_20221202T123056.dat !
%    ...
% however, it is very rare to find something with #L(H) < 512

    addpath ../matrix_tool
    addpath ../CHM
    N = 11;

    while 1
	% q = 11
        % R = (randi(q, N, N)-1);
        % X = exp(2j*pi*R/q); % random seed (in hope to find BH(11, q) with q > 11
        % printf("------------------------------------------------\nwait: \"sinkhorn\" tries to converge...\n");
	X = randn(N) + 1j *randn(N);
	Y = sinkhorn(X, kMax);

        if nh(Y)<1e-12 && n1(Y)<1e-12
            isBH(Y, 1000);
            Y = orderCore(dephase(Y));
            d = ud(Y, "S", 1e-8);
            L = size(getUniqueLIMITED(haagerup(Y), 512),1); % filter Y : #L(Y) < 1000
            % printf("#L = %d", L);
            if L<512
                fileName = strcat("Y", int2str(N), "_", int2str(d), "_", int2str(L), "_", datestr(now(), 30), ".dat");
                save(fileName, "*");
                printf("\nsaved as %s\n", fileName);
            else
                printf(".");
            end
        else
            printf("x");
        end
    end % 1

end
