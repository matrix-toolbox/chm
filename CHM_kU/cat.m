function U = cat(d, a, b, c)
% ------------------------------------------------------------------------------
% 2022-06-28 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Quantized cat map.
%
%    (3, 1,2,-2) --> 2-unitary U in BH(9, 18) : d(U) = 4, #L(U) = 9
%    (3, 2,4,-4) --> 2-unitary U in BH(9, 9)  : d(U) = 4, #L(U) = 9
%    (3, 4,8,-8) --> 2-unitary U in BH(9, 9)  : d(U) = 4, #L(U) = 9
%
%    in general one recovers BH(9, 9) and BH(9, 18) with d and #L as above (respectively)
%    for the following input parameters:
%    (3, k, 2*k, -2*k) for k=1,2,4,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31, ... probably infinite sequence of integers with gaps...
%                          gaps: without k=3*m
%    this triplet above probably works for any odd dimension d = 3, 5, 7, ...
%
%    other triplets for d=3:
%    (-8, -7, -8) --> BH(9, 18), defect=4, L=9
%    (-4, -8, -4) --> BH(9, 9),  defect=4, L=9
%    (-5, -1, -4) --> BH(9, 18), defect=4, L=9
%    (-1, -2, -2) --> BH(9, 18), defect=4, L=9
%    (-7, -8, -8) --> BH(9, 18), defect=4, L=9
%    :
%    :
% ------------------------------------------------------------------------------

    N = d^2;
    U = zeros(N, N);

    for j=1:N
    for k=1:N
        U(j, k) = (1/1)*exp((1j*pi/N)*(a*(j^2)+b*(k^2)+c*j*k));
    end, end

end

