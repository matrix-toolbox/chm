function tests = fourier_test
% ------------------------------------------------------------------------------
% 2017-01-17 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Standard and Hermirian Fourier matrix -- "unit test";
%
% >> result = runtests('fourier_test')
% >> rt = table(result)
% ------------------------------------------------------------------------------

    tests = functiontests(localfunctions);

end

% ------------------------------------------------------------------------------

function TOLERANCE = TOLERANCE % FIX this ...
    TOLERANCE = 1e-10;
end

% ------------------------------------------------------------------------------

function testUnitarityCondition(testCase)
    disp('Wait... Checking unitarity condition for ''classic'' Fouriers...')
    MAXIMUM_N = 18;
    for N = 1 : MAXIMUM_N % N = 18 is an upper bound for MATLAB capabilities...
        disp(sprintf('N (of %d) = %d', MAXIMUM_N, N));
        H = fourier(N, 'hermitian');
        T = N * N;
        assert(norm(H * H' - T * eye(T), 'fro') <= TOLERANCE, ...
        'Unitarity condition failed!')
    end
end

% ------------------------------------------------------------------------------

function testHermiticityCondition(testCase)
    disp('Wait... Checking hermiticity condition for ''hermitian'' Fouriers...');
    MAXIMUM_N = 32;
    for N = 1 : MAXIMUM_N
        disp(sprintf('N (of %d) = %d', MAXIMUM_N, N));
        H = fourier(N, 'hermitian');
        assert(sum(abs(imag(eig(H)))) <= TOLERANCE, 'Spectrum is not real!')
        assert(norm(H - H', 'fro') <= TOLERANCE, 'Matrix is not Hermitian!')
    end
end

% ------------------------------------------------------------------------------

function testEquivalenceRelation(testCase)
    disp('Wait... Checking equivalence relation between ''classic'' and ''hermitian'' Fouriers...');
    MAXIMUM_N = 32;
    for N = 1 : MAXIMUM_N
        disp(sprintf('N (of %d) = %d', MAXIMUM_N, N));
        H = fourier(N, 'hermitian');
        F = fourier(N, 'classic');
        P = H' * kron(F, F) / N / N;
        for row = 1 : N * N;
            assert(abs(sum(P(row, :)) - 1) <= TOLERANCE, ...
            'Matrix P is not a permutation matrix: F_n (x) F_n = H_N * P');
        end
    end
end
