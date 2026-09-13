function mfprintTest
% 20171219
% testing routine for mfprint.m
%
% exemplary output for F8_5_11111:
%
%    v_i(2)  m_i(2)
%    0       336
%    2       448
%
%
%    v_i(3)  m_i(3)
%    0       1344
%    4       1792
%
%
%    v_i(4)  m_i(4)
%    0       1428
%    8       3136
%    16      336

  F8_5_11111 = [ % standard Fourier in a non-standard form, cf. arXiv:1204.5160v1
    1  1  1  1  1  1  1  1;
    1  1  1  1 -1 -1 -1 -1;
    1  1 -1 -1  1  1 -1 -1;
    1  1 -1 -1 -1 -1  1  1;
    1 -1  1 -1  1 -1  1 -1;
    1 -1  1 -1 -1  1 -1  1;
    1 -1 -1  1  1 -1 -1  1;
    1 -1 -1  1 -1  1  1 -1;
  ];
  PHI = mfprint(F8_5_11111);

  w = exp(2 * pi * i / 3);
  P7 = [ % Petrescu matrix, cf. PhD Thesis UCLA (1997)
    -1  1  w  1  w  1  w;
     1 -1  w  1  1  w  w;
     w  w -w  1  w  w  1;
     1  1  1 -1  w  w  w;
     w  1  w  w -w  w  1;
     1  w  w  w  w -w  1;
     w  w  1  w  1  1 -1;
  ]; 
  PHI = mfprint(P7);

end

