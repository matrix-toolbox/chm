function H = GH7_Z6
% 2025-11-03
%
% Matrix described by Bradley Brock in 1988 https://www.sciencedirect.com/science/article/pii/0097316588900544
% Later this matrix was used by M. Petrescu who introduced a 1-parameter affine family...
% This matrix also appears in the note: https://arxiv.org/pdf/1510.06816

    w = exp(2j*pi/3);

    H = [
        -1     1   1   1      1   1   1  ;

         1    -w   w   w      w^2 1   1  ;
         1     w  -w   w      1   w^2 1  ;
         1     w   w  -w      1   1   w^2;

         1     w^2 1   1     -w   w   w  ;
         1     1   w^2 1      w  -w   w  ;
         1     1   1   w^2    w   w  -w  ;
    ]

end
