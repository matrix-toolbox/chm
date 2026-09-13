function summary(X, Epsilon=1e-8, get_defect = 0, get_lambda = 0)
% ------------------------------------------------------------------------------
% 2022-08-27
% 2025-05-25
%
% summarize CHM
%
% DEFAULT USAGE: summary(X);
% DETAILED USAGE: summary(CHM, Epsilon, get-defect = 1, get-Lambda = 1);
% ------------------------------------------------------------------------------

    addpath ../matrix_tool

    printf("-----------------------------------------------\n");
    printf("dimensions: [%d %d]\n", size(X, 1), size(X, 2));
    printf("CHM norm = %g\n", nh(X));
    printf("|1| norm = %d\n", n1(X));
    if nh(X)>1e-6 || n1(X) >1e-6, printf("\n\t********\n\t* fail *\n\t********\n\n"); return, end

    if get_defect
        printf("defect   = %d\n", ud(X, "S", 1e-8));
    else
        printf("skipping defect...\n");
    end

    if get_lambda
        printf("#L       = %d\n", cL(X, Epsilon));
    else
        printf("skipping Lambda-set...\n");
    end

    printf("symmetry = %s\n", cSH(X));
    [q verdict] = isBH(X, 10000, "VERBOSE");
%    printf("BH(N, q) = %s\n", verdict);
    printf("-----------------------------------------------\n");

end
