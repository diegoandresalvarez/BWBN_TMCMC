function p_jp1 = calculate_pjp1(log_fD_theta_j, p_j)
%% p_jp1 = calculate_pjp1(log_fD_theta_j, p_j)
%
%   Function used to compute the tempering parameter p_(j+1) in TMCMC
%
%   Find p_(j+1) such that Coefficient of variation of fD_theta_j <= to a
%   prescribed threshold, i.e., define w_j = fD_theta_j.^(pj1 - pj), then:
%
%                                 std(w_j)
%             COV(fD_theta_j) = ------------ <= threshold 
%                                 mean(w_j)
%
%   (*) Input data:
%
%   - log_fD_theta_j: Log. of fD_theta_j
%   - p_j:            Value of p_(j)
%
%   (*) Output data:
%
%   - p_jp1: Value of p_(j+1)
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            04-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning
threshold = 1;          % 100% = threshold for the Cefficient of variation

% w_j = @(e) fD_theta_j^e;                % N x 1 vector
% Note the following trick in order to compute e = p_(j+1) - p_(j):
% Take into account that e >= 0
w_j = @(e) exp(abs(e)*log_fD_theta_j);     % N x 1 vector

% fmin = @(e) std(w_j(e))/mean(w_j(e)) - threshold;
fmin = @(e) std(w_j(e)) - threshold*mean(w_j(e)) + eps;

% 'fzero' let me find a zero of a function
e = abs(fzero(fmin, 0));      % e is >= 0, and fmin is an even function

if isnan(e)
  % log(realmax) = 709 \approx 500
  % if (p_jp1 - p_j)*max(log_fD_theta_j) > log(realmax)
  % 	the algorithm diverges
  % end
  % therefore:
  % (p_jp1 - p_j)*max(log_fD_theta_j) = 500
  % and 
  % p_jp1 = p_j + 500/max(log_fD_theta_j)
  % that is:
  p_jp1 = min(1, p_j + 500/max(log_fD_theta_j));
  fprintf('Variable p was set to %f, since it is not possible to find a suitable value of p\n',p_jp1);
else
  p_jp1 = min(1, p_j + e);
end

%{
figure
p = linspace(0,0.3,10000); 
hold on
plot(p, arrayfun(fmin, p));
plot(e,0,'rx');
grid minor;
%}

end