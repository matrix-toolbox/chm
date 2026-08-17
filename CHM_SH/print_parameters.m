function print_parameters(p)
% ------------------------------------------------------------------------------
% 2023-01-15 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Module for "solver.m".
% ------------------------------------------------------------------------------

    printf("parameters:\n");
    for j = 1 : size(p, 2)
        rp = real(p(j));
        ip = imag(p(j));
        if rp<0 && ip<0
            printf("p%d = -%2.15g - %2.15gi\n", j, -rp, -ip);
        elseif rp<0 && ip>0
            printf("p%d = -%2.15g + %2.15gi\n", j, -rp, ip);
        elseif rp>0 && ip<0
            printf("p%d = +%2.15g - %2.15gi\n", j, rp, -ip);
        else
            printf("p%d = +%2.15g + %2.15gi\n", j, rp, ip);
        end
    end

end