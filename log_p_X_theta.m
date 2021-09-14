function MSE = log_p_X_theta(theta, X, u, w0, tt, dt, env)
%% MSE = log_p_X_theta(theta, X, u, w0, tt, dt, env)
%
%  This function computes the log-likelihood PDF of the Mean Square Error
%  (MSE) between 'X' and the estimated response of the system given the
%  parameters 'theta'.
%
%  BIBLIOGRAPHY:
%
% - MUTO, Matthew; BECK, James L. "Bayesian updating and model class
%   selection for hysteretic structural models using stochastic
%   simulation". Journal of vibration and control. 14(1-2):7-34, 2008.
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            11-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning
N   = length(u);
N0  = size(theta,1);
MSE = zeros(N0,1);         % Allocate space in memory for MSE

for i = 1:N0
  %% Definition of the state function
  BW_test = @(x,u) diff_eq_real(x,u,[w0 theta(i,3:end)]);
  F_BWRK  = @(x_k,u_k) rk_discrete(BW_test,x_k,u_k,dt);

  %% Initial condition:
  % Initial condition x, xd, z, e (Displacement, velocity. hysteretic
  % displacement, hysteretic energy) set to zero
  x_0 = zeros(4,1);

  % vector where I'm going to store the system response
  x_kk      = zeros(4,N);
  x_kk(:,1) = x_0;                   % Initial state

  %% Computing system response
  for j = 2:N
    x_kk(:,j) = F_BWRK(x_kk(:,j-1),u(j-1));
  end

  x_kk = x_kk';

  %% Total dissipated energy
  % Compute 'dissipated elastic energy'
  %
  %                                            /t_f
  %  dissipated_elastic_energy = alpha * w0^2 *|     est_displ * est_vel dt
  %                                            /t_0
  %
  % Rememeber that this is a cummulative measure. The eq. is multiplied by
  % 1e-6 because of the units, without the factor, the eq. has units of
  % 'mm^2/s^2'; with the factor, the eq. (dissipated hysteretic energy) has
  % units of 'J/kg'.
  %
  alpha = theta(i,4);
  diss_elastic_energy = 1e-6*alpha*(w0^2)*cumtrapz(tt, x_kk(:,1).*x_kk(:,2));

  % tot_diss_energy = diss_elastic_energy + diss_hysteretic_energy
  tot_diss_energy = diss_elastic_energy + x_kk(:,4);

  %% Compute Mean Square Error (MSE) (MUTO, eq. 19) - Modified
  % Prediction error displacements, divided by envelope for normalization
  % between [-1,1]
  p_error_dis = (x_kk(:,1) - X(:,1))./env(:,1);

  % This step implicitly guarantees that inequalities typically used to
  % ensure the Bounded Input-Bounded Output stability of the
  % BWBN-differential equation are satisfied
  % (See eq. in line 72 of 'diff_eq_real.m' file)
  p_error_dis(imag(p_error_dis) ~= 0) = realmax;

  % standard deviation noise (displacement)
  sigma_q_dis = sqrt(theta(i,1));

  % Prediction error dissipated energy, divided by envelope for
  % normalization between [-1,1]
  p_error_ene = (tot_diss_energy - X(:,2))./env(:,2);

  % This step implicitly guarantees that inequalities typically used to
  % ensure the Bounded Input-Bounded Output stability of the
  % BWBN-differential equation are satisfied
  % (See eq. in line 72 of 'diff_eq_real.m' file)
  p_error_ene(imag(p_error_ene) ~= 0) = realmax;

  % standard deviation noise (energy)
  sigma_q_ene = sqrt(theta(i,2));

  % I add 'realmin' in order to avoid 'Inf' when the PDF = 0:
  %
  %           log(0) = Inf --------> log(0 + realmin) ~= Inf
  %
  MSE(i)  = sum(log(normpdf(p_error_dis, 0, sigma_q_dis) + realmin)) + ...
            sum(log(normpdf(p_error_ene, 0, sigma_q_ene) + realmin));
end

end
%% END
