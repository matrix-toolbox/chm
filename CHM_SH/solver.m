function Y = solver(MetaData)
% ------------------------------------------------------------------------------
% 2023-01-22 Wojciech Bruzda; name[at]uj.edu.pl : name = w.bruzda
%            https://chaos.if.uj.edu.pl/~karol/hadamard/
%            https://github.com/matrix-toolbox/
% ------------------------------------------------------------------------------
% Input data:
%    MetaData{1} = max. number of iterations
%    MetaData{2} = number of independent phases
%    MetaData{3} = tolerance for being CHM
%    MetaData{4} = matrix pattern
%    MetaData{5} = scaling factor (~random walk scope)
% ------------------------------------------------------------------------------
% Example of call:
% octave:181> Y = solver(nthargout([1:5], @HH_8_5_102)); summary(Y, 1e-8);
% ..............................................................................
% Example of output for HH_8_5_102:
% octave:182> Y = solver(nthargout([1:5], @HH_8_5_102)); summary(Y, 1e-8);
% RESET
% RESET
% RESET
% RESET
% RESET
% RESET
% RESET
% RESET
% RESET
% parameters:
% p1 = -0.28030447487089 + 0.95991114243317i
% p2 = -0.90942583111813 + 0.41586615358202i
% p3 = -0.02369649067919 - 0.99971919874007i
% p4 = -0.84395409130978 - 0.53641540969801i
% p5 = +0.96680671909618 - 0.25550884115911i
% p6 = -0.18445964867304 - 0.98284008771082i
% p7 = -0.83054080097954 + 0.55695778826431i
% p8 = +0.01500044965257 + 0.99988748692551i
% p9 = +0.81885166223632 - 0.57400518747901i
% -----------------------------------------------
% CHM norm = 3.79854e-13
% |1| norm = 0
% defect   = 5
% #L       = 102
% symmetry = matrix is Hermitian: X = X^*
% BH(N, q) = probably NOT BH(N=8, q) for 1 < q < 10000...
% -----------------------------------------------
% ------------------------------------------------------------------------------

    addpath ../CHM

    iiMax = MetaData{1};
    nP = MetaData{2}; % number of phases
    do
        printf("RESET\n");
        Z = +Inf;
        Z_OPTIMAL = +Inf;
        mu = 1;
        PHASE = rand(1, nP);
        ii = 0;
        while Z_OPTIMAL > MetaData{3} && ii<iiMax
            ii++;
            RESTORE_PHASE = PHASE;
            PHASE = mod(PHASE + randn(1, nP)*mu, 1);
            for j = 1 : nP
                p(j) = exp(2j*pi*PHASE(j));
            end

            % ------------------------------------------------------------------
            % here put fixed phases, if any...
            % ------------------------------------------------------------------

            Y = MetaData{4}(p);
            Z = nh(Y);
            if Z < Z_OPTIMAL
                Z_OPTIMAL = Z;
                mu = Z * MetaData{5};
                printf("%2.14g \t k=%d\n", Z, ii);
                ii = 0;
            else
                PHASE = RESTORE_PHASE;
            end
        end % while Z...
    until (ii!=iiMax)
    print_parameters(p);
end
