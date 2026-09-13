function nuf=sfc(Y, d_mode="VERBOSE")
% 2022-08-09
% for every matrix element y find all others of the form: y, -y, 1/y, -1/y, i*y, -i*y, i/y, -i/y, y/i, -y/i
% sfc = simple fraction correlations
%
% 2022-08-20
% returns the number of "unknown" fields (for which no simple dependency was detected)
%
% OLD NAME:
% >> n = find_siblings(Y)
% >< n = find_siblings(Y, "VERBOSE")


    d = size(Y, 1);
    Y = Y(2:d,2:d); % get core
    ET = 1e-12;
    YL = reshape(transpose(Y), 1, (d-1)^2); % row by row NOT column by column!!


    excluded_j=[];
    jj = [];
    % find repetitions
    for j=1:(d-1)^2-1
        if !any(excluded_j(:) == j) % check only allowed indices
            y = YL(j);
            for k=j+1:(d-1)^2
                if abs(YL(k)-y)<ET
                    excluded_j = [excluded_j j];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = Y(%d, %d) \t\t= %g + i*%g", j, x1, y1, x2, y2, real(YL(j)),imag(YL(j))));
                    continue
                end
                if abs(YL(k)+y)<ET
                    excluded_j = [excluded_j k];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = -Y(%d, %d) \t\t= %g + i*%g", j, x1, y1, x2, y2, real(YL(j)),imag(YL(j))));
                    continue
                end
                if abs(YL(k)-1/y)<ET
                    excluded_j = [excluded_j k];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = Y*(%d, %d) \t\t= %g + i*%g", j, x1, y1, x2, y2, real(YL(j)),imag(YL(j))));
                    continue
                end
                if abs(YL(k)+1/y)<ET
                    excluded_j = [excluded_j k];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = -Y*(%d, %d) \t\t= %g + i*%g", j, x1, y1, x2, y2, real(YL(j)),imag(YL(j))));
                    continue
                end
            end
        end
    end

    symbol = ["abcdefghjkmnopqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ"]; % excluded: i, j, I
%    symbol = ["abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"]; % full Latin set
    map_j_symbol = {};
    for kk=1:size(unique(jj),2)
%        disp(sprintf("%d --> %c", sort(unique(jj))(kk), symbol(kk)));
        map_j_symbol{sort(unique(jj))(kk)} = symbol(kk);
    end
%   now we have, eg.: map_j_symbol{3} returns "b"


%   2nd iteraton: create array of symbols
    for xj=1:d
    for yj=1:d
        array_of_symbols{xj}{yj}=".";
        if (xj==1) || (yj==1)
            array_of_symbols{xj}{yj}="1";
        end
    end, end

    excluded_j=[];
    jj = [];
    % find repetitions
    for j=1:(d-1)^2-1
        if !any(excluded_j(:) == j) % check only allowed indices
            y = YL(j);
            for k=j+1:(d-1)^2
                if abs(YL(k)-y)<ET
                    excluded_j = [excluded_j j];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = %c \t Y(%d, %d) = %c", j, x1, y1, map_j_symbol{j}, x2, y2, map_j_symbol{j}));
                    array_of_symbols{x1}{y1}=map_j_symbol{j};
                    array_of_symbols{x2}{y2}=map_j_symbol{j};
                    continue
                end
                if abs(YL(k)+y)<ET
                    excluded_j = [excluded_j k];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = %c \t Y(%d, %d) = -%c", j, x1, y1, map_j_symbol{j}, x2, y2, map_j_symbol{j}));
                    array_of_symbols{x1}{y1}=map_j_symbol{j};
                    array_of_symbols{x2}{y2}=["-", map_j_symbol{j}];
                    continue
                end
                if abs(YL(k)-1/y)<ET
                    excluded_j = [excluded_j k];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = Y*(%d, %d) \t\t= %g + i*%g", j, x1, y1, x2, y2, real(YL(j)),imag(YL(j))));
                    array_of_symbols{x1}{y1}=map_j_symbol{j};
                    array_of_symbols{x2}{y2}=["1/", map_j_symbol{j}];
                    continue
                end
                if abs(YL(k)+1/y)<ET
                    excluded_j = [excluded_j k];
                    excluded_j = [excluded_j k];
                    jj = [jj j];
                    [x1 y1] = jc2jf(j, d);
                    [x2 y2] = jc2jf(k, d);
%                    disp(sprintf("%d\tY(%d, %d) = -Y*(%d, %d) \t\t= %g + i*%g", j, x1, y1, x2, y2, real(YL(j)),imag(YL(j))));
                    array_of_symbols{x1}{y1}=map_j_symbol{j};
                    array_of_symbols{x2}{y2}=["-1/", map_j_symbol{j}];
                    continue
                end
            end
        end
    end


    % 3rd iteration: find +1, -1, +1j, -1j[, and simple omegas...]
    for q=3:12
    for p=3:12
        omega = exp(2j*pi/q);
        for xj=1:(d-1)
            for yj=1:(d-1)
                as = Y(xj, yj);
                if abs(as-1.0)<ET, array_of_symbols{xj+1}{yj+1}="1"; end
                if abs(as+1.0)<ET, array_of_symbols{xj+1}{yj+1}="-1"; end
                if abs(as-1j)<ET, array_of_symbols{xj+1}{yj+1}="1j"; end
                if abs(as+1j)<ET, array_of_symbols{xj+1}{yj+1}="-1j"; end
                %if abs(as-omega^p)<ET, array_of_symbols{xj+1}{yj+1}="omega"; end
            end
        end
    end, end





#    % 4th iteration: find 1j-functions: Y(j)/Y(j+Delta_j) = +/-1j
#    for j=1:(d-1)^2-1
#        y = YL(j);
#        for k=j+1:(d-1)^2
#            y2 = YL(k);
#            if abs(y/y2-1j)<ET
#                 [xj yj] = jc2jf(j, d);
#                 [xk yk] = jc2jf(k, d);
#                 array_of_symbols{xj}{yj} = [array_of_symbols{xk}{yk}, "*1j"];
#            end
#        end
#    end


    if (!strcmp(d_mode, "SILENT"))
        % display symbolic matrix
        printf("\n\n-------------------------------------------------------------------------------");
        printf("\n");
        for xj=1:(d)
            printf("\n");
            for yj=1:(d)
                printf("%s\t", array_of_symbols{xj}{yj});
            end
        end
        printf("\n\n-------------------------------------------------------------------------------");
        printf("\n\n");
    endif

    if (strcmp(d_mode, "VERBOSE"))
        % 5th iteration: find (additional) simple (non-trivial) functional dependences
        for j=1:(d-1)^2-2
            y = YL(j);
            for k=j+1:(d-1)^2-1
                y2 = YL(k);
                if abs(y2-1.0)<ET || abs(y2+1.0)<ET || abs(y/y2-1.0)<ET || abs(y/y2+1.0)<ET
                    continue % skip trivial fractions
                end
                for l=k+1:(d-1)^2
                    if abs(y/y2-YL(l))<ET
                        [x1 y1] = jc2jf(j, d);
                        [x2 y2] = jc2jf(k, d);
                        [x3 y3] = jc2jf(l, d);
    #                    disp(sprintf("%d\tY(%d, %d)/Y(%d, %d) = Y(%d, %d) \t\t= %g + i*%g", j, x1, y1, x2, y2, x3, y3, real(YL(l)),imag(YL(l))));
                        if array_of_symbols{x1}{y1} == ".", s1 = ["Y(", int2str(x1),",",int2str(y1),")"]; else s1 = array_of_symbols{x1}{y1}; end;
                        if array_of_symbols{x2}{y2} == ".", s2 = ["Y(", int2str(x2),",",int2str(y2),")"]; else s2 = array_of_symbols{x2}{y2}; end;
                        if array_of_symbols{x3}{y3} == ".", s3 = ["Y(", int2str(x3),",",int2str(y3),")"]; else s3 = array_of_symbols{x3}{y3}; end;
                        disp(sprintf("%s/%s = %s \t\t= %g + i*%g", s1, s2, s3, real(YL(l)),imag(YL(l))));
                        continue
                    end
                end
            end
        end
    endif

    nuf = 0;
    for xj=1:(d)
        for yj=1:(d)
            if array_of_symbols{xj}{yj} == ".", nuf++; end
        end
    end


end




function [x, y]=jc2jf(j, d)
% recalculate linear core-index to 2-d full-matrix-index
    x = 1 + 1 + floor((j - 1) / (d - 1));
    y = 1 + 1 + mod(j - 1, d - 1);
end

