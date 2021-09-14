function theta_tp1 = Metropolis_Hastings(theta_t,p,q,sampling_proposal_PDF,burnin)
%% theta_tp1 = Metropolis_Hastings(theta_t,p,q,sampling_proposal_PDF,burnin)
%
% Metropolis-Hastings algorithm.
%
% * Input data:
%
%   - theta_t :              Current value of theta.
%   - p:                     We want to draw samples from the PDF p(theta)
%   - q:                     Proposal PDF
%   - sampling_proposal_PDF: Candidate generating distribution
%   - burnin:                Burn-in (optional)
%
% * Output data:
%
%   - theta_tp1: accepted sample
%
% BIBLIOGRAPHY:
%
% - WALSH, B. "Markov Chain Monte Carlo and Gibbs sampling". Lecture notes
%   for EEB 581. April 26, 2004.
% - ANDRIEU, Christophe, et. al. "An introduction to MCM for machine
%   learning". Machine learning, 50, 5-43, 2003.
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            13-Aug-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning

%% In case of no burn-in and/or thinning (Lag)
if nargin < 5
  burnin = 1;               % No burn-in
end

%% Burn-in
for i = 1:burnin            % Carry out burn-in
  % Sampling from candidate PDF
  theta_ast = sampling_proposal_PDF(theta_t);

  % Compute alpha = Ratio of the density at the candidate (theta_ast) and
  % current (theta_t) points
  alpha = (p(theta_ast)*q(theta_t,theta_ast))/...
                                (p(theta_t)    *q(theta_ast,theta_t));

  % Metropolis-Hastings criterion
  if rand <= min(alpha,1)
    theta_tp1 = theta_ast;   % Accept the candidate
  else
    theta_tp1 = theta_t;     % Reject the candidate and use the same state
  end
end

end
%% END