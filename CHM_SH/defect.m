function defect(H)
% 2023-02-18

    N = size(H, 1);
    m = (N + 1)*N/2; % #variables = #columns in R
    c = (N - 1)*N/2; % #row in R

    R = zeros(c, m);
%    size(R)

    r = 0; % row of R
    for a=1:N % row of H
    for b=a:N % column of H
        if a==b, continue, end % exclude diagonal
        r++;
        printf("r = %d\n", r);
        for c=1:N % inner product
            j1 = a;
            j2 = c;
            j3 = c;
            j4 = b;
            if j1>j2; t=j1; j1=j2; j2=t; end;
            if j3>j4; t=j3; j3=j4; j4=t; end;
            printf(" H(%d, %d)*conj(H(%d, %d))", j1, j2, j3, j4);


            offset = j2-(j1-1);
            for m=1:j1-1
                offset += (N-m+1);
            end
            printf(" R(%d, %d)", r, offset);
            R(r, offset) += H(j1, j2)*conj(H(j3, j4));


            offset = j4-(j3-1);
            for m=1:j3-1
                offset += (N-m+1);
            end
            printf(" -R(%d, %d)", r, offset);
            R(r, offset) += -H(j1, j2)*conj(H(j3, j4));
           

            printf("\n");
        end        
    end, end
    R = [real(R); imag(R)];
    size(R)

    
    m = (N + 1)*N/2; % #variables = #columns in R
    r = sum(svd(R) > 1e-8); % rank(R);
    p = N;
    printf("#variables = %d\n", m);
    printf("   rank(R) = %d\n", r);
    printf(" dephasing = %d\n", p);
    printf("    defect = %d\n", m - r - p);

    UNITARY_DEFECT = ud(H, "S", 1e-8)   
    
end

