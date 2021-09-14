function [theta_m,log_S,p,Theta] = tmcmc(log_fD_theta,f_theta,sample_from_f_theta,N,burnin,last_burnin)
%% [theta_m,log_S,p,Theta] = tmcmc(log_fD_theta,f_theta,sample_from_f_theta,N)
%
% Transitional Markov Chain Monte Carlo algorithm.
%
% (*) Input data:
%
%     - log_fD_theta:        Log. of Likelihood PDF (handle function)
%     - f_theta:             Prior PDF (handle function)
%     - sample_from_f_theta: Function that samples from the prior PDF
%                            f_theta (handle function)
%     - N:                   Number of samples to generate each stage. The
%                            algorithm assumes that N0 = N1 = ... = Nm.
%
%                            Parameters of the Metropolis-Hastings algorithm
%     - burnin:
%     - last_burnin:         burn-in in the last iteration
%
% (*) Output data:
%
%     - theta_m: Estimated parameters (N x D matrix)
%     - log_S:   Estimated log. of evidence
%     - p:       Tempering parameters
%     - Theta:   Evolution of parameters (cell array)
%
% BIBLIOGRAPHY:
%
% - CHING, Jianye; CHEN, Yi-Chun. "Transitional Markov Chain Monte Carlo
%   method for Bayesian model updating, model class selection and model
%   averaging". Journal of Engineering Mechanics. ASCE. 133(7):816-832,
%   July 1, 2007.
%
% - MUTO, Matthew; BECK, James L. "Bayesian updating and model class
%   selection for hysteretic structural models using stochastic
%   simulation". Journal of vibration and control. 14(1-2):7-34, 2008.
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            05-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Plot graphics while TMCMC is running (better when number of parameters
%  to be estimated = 2)
PLOT_GRAPHICS = false;      % 'true' = plot; 'false' = do not plot

%% I use this line during resampling stage
WITH_REPLACEMENT = true;    % DO NOT CHANGE!!!

%% Constants
% beta is a prescribed scaling factor that suppress the rejection rate and,
% at the same time, makes large MCMC jumps (CHING, eq. 17).
beta   = 0.2;
est_it = 500;                     % Estimated number of iterations
S      = ones(est_it,1);          % Allocate space for S_j (CHING, eq. 15)

%% Obtain N samples from the prior pdf f(T)
j       = 0;                      % Initialize counter
theta_j = sample_from_f_theta(N); % theta_0 (N x D matrix)
p_j     = 0;                      % p_0 = 0 (initial tempering parameter)
D       = size(theta_j, 2);       % size of vector theta

%% Allocate space in memory for theta_jp1
theta_jp1 = zeros(N, D);

%% Cell where I'm going to store all estimated theta matrices during TMCMC
Theta    = cell(1,est_it);
Theta{1} = theta_j;               % Store initial state

%% Vector where I'm going to store tempering parameters
p = zeros(est_it,1);

%% NOTE:
%  I defined a cell array of 'est_it' elements for 'Theta', and a vector
%  with 'est_it' elements for 'p'. If the algorithm needs more than
%  'est_it' iterations, you must increase that number in order to store all
%  the data for posterior Analysis. At the end of the algorithm, the
%  unnecesary space is deleted.

%% TMCMC
while p_j < 1
  %% Plot graphics that show the evolution of estimated parameters
  %  (sampled points).
  if PLOT_GRAPHICS
    figure
    subplot(2,2,2);
    hist(theta_j(:,1), ceil(sqrt(N)));          % histogram
    xlabel('\theta_{1}')
    ylabel('Frequency')
    subplot(2,2,4);
    hist(theta_j(:,2), ceil(sqrt(N)));          % histogram
    xlabel('\theta_{2}')
    ylabel('Frequency')
    subplot(2,2,[1 3]);
    plot(theta_j(:,1), theta_j(:,2), 'b.');     % Sampled points
    hold on;
    ax = axis;
    [xx, yy] = meshgrid(linspace(ax(1),ax(2),100),...
                        linspace(ax(3), ax(4), 99));
    if j == 0
       zz = reshape(f_theta([xx(:) yy(:)]), 99, 100);
    else
       zz = reshape(f_jp1([xx(:) yy(:)]), 99, 100);         
    end
    % contour plots for the intermediate (prior) PDF being plotted
    contour(xx, yy, zz, 50, 'r');
    grid on;
    xlabel('\theta_{1}')
    ylabel('\theta_{2}')
    title(sprintf(...
       'Samples of f_{%d} and contour levels of f_{%d} (red) and f_{%d} (black)', ...
       j, j, j+1));
  end

  j = j+1;                          % Increase counter

  %% Compute tempering parameter p_(j+1)
  log_fD_theta_j = log_fD_theta(theta_j);
  p_jp1          = calculate_pjp1(log_fD_theta_j, p_j);
  p(j+1)         = p_jp1;

  % Print information in the screen while TMCMC is running
  fprintf('TMCMC: Iteration j = %3d, p_(j+1) = %f\n', j, p_jp1);
  
  %% Compute plausibility weights w(theta_j) (CHING, eq. 12) and factors
  %  'S_j' for the evidence (CHING, eq. 15)
  w_j      = exp((p_jp1 - p_j)*log_fD_theta_j);
  w_j_norm = w_j./sum(w_j);         % normalized weights 
  S(j)     = mean(w_j);

  %% Resampling step to obtain N samples from f_{j+1}(theta). Then perform
  %  Metropolis-Hastings on each of these samples using as a 
  %  stationary PDF "f_{j+1}(theta)":
  %  f_jp1 = @(t) f_theta(t).*fD_theta(t).^(p_jp1);   % CHING, eq. 11
  log_f_jp1 = @(t) log(f_theta(t)) + p_jp1*log_fD_theta(t);

  if PLOT_GRAPHICS
    % In the definition of f_jp1 we are including the normalization
    % constant prod(S(1:j)).
    f_jp1 = @(t) exp(log(f_theta(t)) + p_jp1*log_fD_theta(t) - ...
                                                    sum(log(S(1:j))));
    zz = reshape(f_jp1([xx(:) yy(:)]), 99, 100);
    contour(xx, yy, zz, 50, 'k');
  end

  % Use as proposal PDF a Gaussian centered at theta_j(idx,:) and 
  % with covariance matrix equal to an scaled version of the covariance
  % matrix of f_{j+1}(theta):   
   
  % weighted mean
  mu = zeros(1, D);
  for l = 1:N
    mu = mu + w_j_norm(l)*theta_j(l,:);      % 1 x N matrix
  end

  % Compute scaled covariance matrix of f_{j+1}(theta) (CHING, eq. 17)
  Sigma_j = zeros(D);
  for k = 1:N
    % I had to change slightly this formula due to the size of the vectors,
    % and because Ching and Chen forgot to normalize the weight w_j:
    tk_mu   = theta_j(k,:) - mu;
    Sigma_j = Sigma_j + w_j_norm(k)*(tk_mu' * tk_mu);   % CHING, eq. 17
  end
  Sigma_j = beta^2 * Sigma_j;

  % Proposal distribution
  prop_pdf = @(x,y) mvnpdf(x, y, Sigma_j); %q(x,y) = q(x|y).
  prop_rnd = @(x) mvnrnd_box(x, Sigma_j, 1, f_theta);
  
  %% During the last iteration we require to do a better burnin in order 
  % to guarantee the quality of the samples:
  if p_jp1 == 1
    burnin = last_burnin;
  end

  %% Start N different Markov chains
  % Parallel FOR-loop
  parfor i = 1:N
    %% Sample one point with probability w_j_norm
    idx = randsample(N, 1, WITH_REPLACEMENT, w_j_norm);

    %% Metropolis-Hastings routine (Matlab's routine)
    % smpl = mhsample(start, nsamples, 'pdf', pdf,
    %                 'proppdf', prop_pdf, 'proprnd', prop_rnd);
    % start    = row vector containing the start value of the Markov Chain,
    % nsamples = number of samples to be generated
    theta_jp1(i,:) = mhsample(theta_j(idx, :), 1, 'logpdf',  log_f_jp1, ...
                                                  'proppdf', prop_pdf, ...
                                                  'proprnd', prop_rnd, ...
                                                  'burnin',  burnin);

    %{
    [theta_jp1(i,:), acceptance_rate] = mhsample(theta_j(idx, :), 1, ...
                                                 'logpdf',  log_f_jp1, ...
                                                 'proppdf', prop_pdf, ...
                                                 'proprnd', prop_rnd, ...
                                                 'burnin',  burnin);
    %}

    % According to Cheung and Beck (2009) - Bayesian model updating ...,
    % the initial samples from reweighting and the resample of samples of
    % f_j, in general, do not exactly follow f_jp1, so that the Markov
    % chains must "burn-in" before samples follow f_jp1, requiring a large
    % amount of samples to be generated for each level.

    %{
    %% Adjust the acceptance rate (optimal = 23%)
    %    See: http://www.dms.umontreal.ca/~bedard/Beyond_234.pdf
    if acceptance_rate < 0.3
      % Many rejections means an inefficient chain (wasted computation
      % time), decrease the variance
      beta = 0.99*beta;
    elseif acceptance_rate > 0.5
      % High acceptance rate: Proposed jumps are very close to current
      % location, increase the variance
      beta = 1.01*beta;
    end
    %}
  end

  %% Prepare for the next iteration
  Theta{j+1} = theta_jp1;
  theta_j    = theta_jp1;
  p_j        = p_jp1;
end

%% TMCMC provides N samples distributed according to f(T|D)
theta_m = theta_j;

%% Compute evidence (normalization constant in Bayes' theorem)
if j < est_it-2
  S(j+2:end) = [];
end
log_S = sum(log(S(1:j)));

%% delete unnecesary data
if j < est_it-2
  p(j+2:end)     = [];
  Theta(j+2:end) = [];
end

end
%% END
