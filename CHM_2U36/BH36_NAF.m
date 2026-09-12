function [H, info] = BH36_NAF(t)
% ------------------------------------------------------------------------------
% 2026-08-24
% analytic non-affine deformation of the two-unitary BH(36, 6)
% ------------------------------------------------------------------------------
% >>  H = BH36_NAF(t);
% >>  [H info] = BH36_NAF(t)
% >>  B = BH36_NAF(0); % compare with A = A36();
%
% it returns the unnormalised complex Hadamard matrix H (so H*H' = 36*I when the ordinary flattening is used)
% -- divide by 6 for the unitary normalisation
%
% the construction is the one-parameter non-affine branch through the Butson matrix in 10.1515/spma-2024-0010
% the essential free parameter is one real phase t with x = exp(i*t)
% it is periodic modulo 2*pi but only the following arc gives unimodular y and z:
%
%   t in [t_L, t_R] modulo 2*pi
%
%   t_L = -4.48651889056550858...
%   t_R =  0.29772868577911760...
%
% equivalently, for the principal representative t in [-pi, pi]:
%
%   t in [-pi, 0.29772868577911760...]
%        union
%        [1.79666641661407790..., pi]
%
% the two endpoints are allowed but the quadratic discriminant vanishes there, so the parametrisation is not differentiable at those points
% values in the complementary open arc are currently rejected
%
% the optional output INFO contains the reduced angle and the three phases x, y, z used in the multiplier
%
%   [H, p] = BH36_NAF(0.1);
%   U = H/6;
%   disp([p.x, p.y, p.z]);
% ------------------------------------------------------------------------------
  if nargin ~= 1
    error('BH36_NAF requires exactly one real parameter t!');
  end
  if ~isnumeric(t) || ~isscalar(t) || ~isreal(t) || ~isfinite(t)
    error('parameter (phase) t must be a finite real scalar!');
  end

  % the two unit-circle zeros of the quadratic discriminant
  theta_plus   = 0.2977286857791175962867078421882024;
  theta_second = 1.7966664166140778960217210799981328;
  t_L = theta_second - 2*pi;
  t_R  = theta_plus;

  % use the representative in [t_L, t_L + 2*pi)
  % the admissible connected branch through t=0 is [t_L, t_R]
  tr = mod(t - t_L, 2*pi) + t_L;
  domain_tol = 5e-13;
  if tr > t_R + domain_tol
    error(['the input angle is outside the unimodular branch, modulo 2*pi, ', ...
           'use t in [%.16g, %.16g], in [-pi,pi], use ', ...
           '[-pi, %.16g] union [%.16g, pi]'], ...
           t_L, t_R, t_R, theta_second);
  end

  % clamp tiny endpoint overshoots caused by floating-point reduction
  if tr > t_R
    tr = t_R;
  end
  if tr < t_L
    tr = t_L;
  end

  zeta = exp(1i*pi/3);
  x = exp(1i*tr);

  % F(x,z) = A(x) z^2 + B(x) z + C(x)
  Acoef = 1 + zeta + (zeta - 2)*x;
  Bcoef = x^2 - 4*zeta*x + zeta^2;
  Ccoef = x*((2*zeta - 1)*x + 2 - zeta);
  Delta = Bcoef^2 - 4*Acoef*Ccoef;

  % select the analytic square-root sheet fixed by sqrt(Delta(1)) = 2-zeta
  % a short phase-unwrapping path from 0 to tr avoids discontinuities of the principal complex square root
  discriminant_tol = 2e-12;
  if abs(Delta) <= discriminant_tol
    z = -Bcoef/(2*Acoef); % repeated root at an endpoint
  else
    ngrid = max(17, ceil(abs(tr)/0.01) + 1);
    tg = linspace(0, tr, ngrid);
    xg = exp(1i*tg);
    Ag = 1 + zeta + (zeta - 2).*xg;
    Bg = xg.^2 - 4*zeta.*xg + zeta^2;
    Cg = xg.*((2*zeta - 1).*xg + 2 - zeta);
    Dg = Bg.^2 - 4.*Ag.*Cg;
    phi = unwrap(angle(Dg));
    sqrtDelta = sqrt(abs(Delta))*exp(0.5*1i*phi(end));
    z = (-Bcoef - sqrtDelta)/(2*Acoef);
  end

  denom = x*(z + 2*zeta - 1) - zeta*z;
  if abs(denom) < 1e-12
    error('the rational formula for y is singular at this input!');
  end
  y = zeta*x*z/denom;

  % guard against using the formula outside its CHM arc or against a branch-selection failure in very low precision
  phase_error = max(abs(abs([x, y, z]) - 1));
  if phase_error > 5e-8
    error(['the computed multiplier is not unimodular (maximum modulus ', ...
           'error %.3e), the input is outside the valid branch or the ', ...
           'working precision is insufficient!'], phase_error);
  end

  B = bh36_exponent_matrix();
  H0 = zeta.^B;

  Phi = [y, 1; ...
         1, z; ...
         1, x; ...
         x, z; ...
         x, y; ...
         y, z];

  Z = zeros(36,36);

  % symbolic labels used only to print the explicit Z(x,y,z) pattern
  Phi_symbol = ['y1'; ...
                '1z'; ...
                '1x'; ...
                'xz'; ...
                'xy'; ...
                'yz'];
  Z_pattern = repmat(' ', 36, 36);

  for aa = 0:5
    for bb = 0:5
      row = 6*aa + bb + 1;
      for cc = 0:5
        s = mod(aa + cc, 6) + 1;
        for dd = 0:5
          col = 6*cc + dd + 1;
          ss = mod(bb + dd, 2) + 1;
          Z(row,col) = Phi(s, ss);
          Z_pattern(row,col) = Phi_symbol(s, ss);
        end
      end
    end
  end

  % print the symbolic 36x36 multiplier pattern in six 6-column blocks
  fprintf('\nexplicit symbolic pattern of Z(x,y,z):\n\n');
  for row = 1:36
    fprintf('%s %s %s %s %s %s\n', ...
            Z_pattern(row,  1: 6), Z_pattern(row,  7:12), ...
            Z_pattern(row, 13:18), Z_pattern(row, 19:24), ...
            Z_pattern(row, 25:30), Z_pattern(row, 31:36));
    if mod(row, 6) == 0 && row < 36
      fprintf('\n');
    end
  end
  fprintf('\n');

  H = H0 .* Z;

  if nargout > 1
    info = struct();
    info.free_parameter_count = 1;
    info.t_input = t;
    info.t_reduced = tr;
    info.x = x;
    info.y = y;
    info.z = z;
    info.discriminant = Delta;
    info.phase_modulus_error = phase_error;
    info.domain_unwrapped = [t_L, t_R];
    info.domain_principal = [-pi, t_R; theta_second, pi];
    info.endpoint_angles = [t_R, theta_second];
    info.normalisation = 'H is unnormalised; use H/6 for a unitary matrix';
    info.branch = 'analytic sheet through x=y=z=1 at t=0';
  end
end


function B = bh36_exponent_matrix()
% exponent matrix B from Eq. (17) / Appendix A of 10.1515/spma-2024-0010
% entries are interpreted modulo 6 and H0 = exp(i*pi*B/3)

  Btxt = [ ...
    '155513424424420002155513040040002420'; ...
    '555131133133044402222404022022135311'; ...
    '335153004004000242335153220220242000'; ...
    '135311313313224042402044202202315551'; ...
    '515333244244240422515333400400422240'; ...
    '315551553553404222042224442442555131'; ...
    '424151531440311402040313113204311402'; ...
    '244511200351555404133400351200222131'; ...
    '220553333242113204442115515000113204'; ...
    '040313002153351200535202153002024533'; ...
    '022355135044515000244511311402515000'; ...
    '442115404555153002331004555404420335'; ...
    '333515224042004004515333224042442442'; ...
    '113531513155202202204440240422313313'; ...
    '351533242000022022533351242000400400'; ...
    '131555531113220220222404204440331331'; ...
    '315551200024040040551315200024424424'; ...
    '155513555131244244240422222404355355'; ...
    '125210010343032123125210454121450305'; ...
    '525434325052432341252101430103341432'; ...
    '305450250523212303305450034301030545'; ...
    '105014505232012521432341010343521012'; ...
    '545030430103452543545030214541210125'; ...
    '345254145412252101012521250523101252'; ...
    '010010143501341105454454501143341105'; ...
    '430430054230525101541541503321252434'; ...
    '412412545303143501250250303545143501'; ...
    '232232450032321503343343305123054230'; ...
    '214214341105545303052052105341545303'; ...
    '034034252434123305145145101525450032'; ...
    '545030254345412145303212254345034301'; ...
    '501410543452010343410501210125505232'; ...
    '503054212303430103321230212303052325'; ...
    '525434501410034301434525234143523250'; ...
    '521012230321454121345254230321010343'; ...
    '543452525434052325452543252101541214'];

  B = double(Btxt) - double('0');
end
