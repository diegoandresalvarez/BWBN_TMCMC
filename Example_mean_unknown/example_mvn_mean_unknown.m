%% Problem statement: (BISHOP, exercise 2.40)
% 
% Consider a D-dimensional gaussian random variable x with PDF
% mvnpdf(x,mu,Sigma), in which the covariance Sigma is known and for which 
% we wish to infer the mean mu from a set of observations:
%
%                    X = {x_1, x_2, ..., x_N}.
%
% Given the prior PDF
%
%                    p(mu) = mvnpdf(x,mu0,Sigma0),
%
% find the corresponding posterior PDF p(mu|X).
%
% Solution:
%
% Likelihood:
%
%               p(X|mu) = prod_{n=1}^N mvnpdf(x_n,mu,Sigma)
%
% Posterior PDF:
% invSigma  = inv(Sigma) 
% invSigma0 = inv(Sigma0)
% maximum likelihood estimator for the mean: mu_ML
% mu_ML      = (1/N)*sum_{n=1}^N x_n
%
%                    p(mu|X) = mvnpdf(x,muN,sigmaN)
%
% where:
% invSigma_N = invSigma0 + N*invSigma;
% Sigma_N    = inv(invSigma_N);
% muN        = Sigma_N*(N*invSigma*mu_ML' + invSigma0*mu0')
%
%
% Here, this problem is addressed analitically and using TMCMC. Results are
% compared in the graphs.
%
% BIBLIOGRAPHY:
%
% - BISHOP, Christopher M. "Pattern recognition and machine learning".
%   Springer. 2006.
% - MURPHY, Kevin P. "Conjugate Bayesian analysis of the Gaussian
%   distribution". Last updataed: October 3, 2007. Retrieved from:
%
%       http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
%
% - CHING, Jianye; CHEN, Yi-Chun. "Transitional Markov Chain Monte Carlo
%   method for Bayesian model updating, model class selection and model
%   averaging". Journal of Engineering Mechanics. ASCE. 133(7):816-832,
%   July 1, 2007.
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            05-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning
close all; clear all; clc;

%% Set Matlab random number stream as 'Mersenne Twister' (mt19937ar)
% s = RandStream('mt19937ar','Seed',0);
% RandStream.setGlobalStream(s);

%% Define the true PDF
true_mu = [ 1 2 ];                % unknown mean (to be estimated)
Sigma   = [ 1    0.4
            0.4  2   ];           % Known covariance matriz

%% Observations
N     = 100;                      % number of observations (x_1,...,x_N)
data  = mvnrnd(true_mu, Sigma, N);
mu_ML = mean(data);               % mean (ML estimate) (BISHOP, Eq. 2.143)

%% Define Prior hyperparameters and prior PDF
mu_0    = [-0.1 -0.6];                  % Mean
Sigma_0 = [ 1    0.5
            0.5  1   ];                 % Variance

prior = @(x) mvnpdf(x, mu_0, Sigma_0);  % Prior PDF

% Function that samples from the prior PDF
p_murnd = @(N) mvnrnd(mu_0, Sigma_0, N);

%% The log-likelihood of X given mu (Sigma is known) - Likelihood PDF
log_p_X_mu = @(mu) ex_log_p_X_mu(mu, data, Sigma);

%% Compute mean and variance of Posterior (analitically)
% Variance (MURPHY, Eq. 211 - BISHOP, Eq. 2.142)
invSigma_0 = inv(Sigma_0);
invSigma   = inv(Sigma);
invSigma_N = invSigma_0 + N*invSigma;
Sigma_N    = inv(invSigma_N);

% Mean (MURPHY, Eq. 212 - BISHOP, Eq. 2.141)
mu_N = Sigma_N*(N*invSigma*mu_ML' + invSigma_0*mu_0');

p_mu_X    = @(x) mvnpdf(x, mu_N', Sigma_N);     % Posterior PDF

% function that samples from the posterior PDF
p_mu_Xrnd = @(N) mvnrnd(mu_N', Sigma_N, N);

%% Plot results obtained analitically

% Figure 1. Input data (observations) + Prior PDF contour plot
figure
subplot(1,2,1);
hold on
plot(data(:,1), data(:,2), 'b.');       % Input data (observations)
ax       = axis;
[xx, yy] = meshgrid(linspace(ax(1),ax(2),100), linspace(ax(3), ax(4), 99));
rr       = reshape(prior([xx(:) yy(:)]), 99, 100);
contour(xx, yy, rr, 50);
grid minor
plot(true_mu(1),true_mu(2),'v', 'MarkerEdgeColor','k',...
                                'MarkerFaceColor','g',...
                                'MarkerSize',10);

plot(mu_ML(1),mu_ML(2),'o', 'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',...
                            'MarkerSize',10);
title('Prior Gaussian PDF and samples from the true PDF');
xlabel('\theta_{1}');
ylabel('\theta_{2}');
legend('Samples of the true PDF (Observations) ~ N(true_mu,Sigma)', ...
       'Prior PDF of \mu', ...
       'Exact value of \mu', ...
       'ML estimate of \mu');

subplot(1,2,2);
hold on
r = p_mu_Xrnd(N);
plot(r(:,1), r(:,2), 'b.');     % Samples from the posterior PDF
ax       = axis;
[xx, yy] = meshgrid(linspace(ax(1),ax(2),100), linspace(ax(3), ax(4), 99));
rr       = reshape(p_mu_X([xx(:) yy(:)]), 99, 100);
contour(xx, yy, rr, 50);
grid minor
axis tight
plot(true_mu(1),true_mu(2),'v', 'MarkerEdgeColor','k',...
                                'MarkerFaceColor','g',...
                                'MarkerSize',10);

plot(mu_ML(1),mu_ML(2),'o', 'MarkerEdgeColor','k',...
                          'MarkerFaceColor','g',...
                          'MarkerSize',10);
title('Posterior Gaussian PDF');
xlabel('\theta_{1}');
ylabel('\theta_{2}');
legend('Samples of the posterior PDF ~ N(mu_N,Sigma_N)', ...
       'Posterior PDF of \mu', ...
       'Exact value of \mu', ...
       'ML estimate of \mu');

%% TMCMC

%  Enable parallel computing within MATLAB
if matlabpool('size') == 0        % Activate all available cores
   matlabpool open
end;

% Start TMCMC. Measure elapsed time.
Nj = 1000;                        % Number of samples per iteration
tic
[theta_fT_D, log_fD, p, Theta] = tmcmc(log_p_X_mu, prior, p_murnd, Nj);
toc

matlabpool close                  % Deactivate all cores

%% Plot TMCMC results

% Figure 1. Tempering parameters
figure;
plot(p);
grid on;
xlabel('Iteration');
ylabel('Tempering parameter');

% Figure 2.
figure
hold on
% Estimated sampled points from the posterior PDF using TMCMC
plot(theta_fT_D(:,1), theta_fT_D(:,2), 'r.');

% True mean
plot(true_mu(1),true_mu(2),'v', 'MarkerEdgeColor','k',...
                                'MarkerFaceColor','g',...
                                'MarkerSize',10);

% maximum likelihood estimator for the mean
plot(mu_ML(1),mu_ML(2),'o', 'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',...
                            'MarkerSize',10);

% Mean of points obtained with TMCMC
mu_TMCMC = mean(theta_fT_D)';
plot(mu_TMCMC(1),mu_TMCMC(2),'^', 'MarkerEdgeColor','k',...
                                  'MarkerFaceColor','g',...
                                  'MarkerSize',10);

ax = axis;
[xx, yy] = meshgrid(linspace(ax(1),ax(2),100), linspace(ax(3), ax(4), 99));
rr       = reshape(p_mu_X([xx(:) yy(:)]), 99, 100);
contour(xx, yy, rr, 50);
grid minor
axis tight
title('Posterior Gaussian PDF + samples from TMCMC');
xlabel('\theta_{1}');
ylabel('\theta_{2}');
legend('Samples of the posterior PDF (provided by TMCMC)', ...
       'Exact value of \mu', ...
       'ML estimate of \mu', ...
       'Mean of points obtained with TMCMC');

%% END