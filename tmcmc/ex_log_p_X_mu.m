function L = ex_log_p_X_mu(mu, X, Sigma)
%% L = ex_log_p_X_mu(mu, X, Sigma)
%
%  This function computes the log-likelihood of a multivariate normal PDF
%  with mean 'mu' and covariance 'Sigma'.
%
%
%% Beginning
n = size(mu,1);
L = zeros(n,1);       % Allocate space in memory for L
for i = 1:n
   L(i) = sum(log(mvnpdf(X, mu(i,:), Sigma)));
end

end
%% END