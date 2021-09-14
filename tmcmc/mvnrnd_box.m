function prop_rnd = mvnrnd_box(mu, sigma, N, box)
%% prop_rnd = mvnrnd_box(mu, sigma, N);
%
%  'N' random vectors from the multivariate normal distribution according
%   to the boundaries defined by 'box'.
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            05-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning

prop_rnd = mvnrnd(mu, sigma, N);

% Keep generating samples until all samples are enclosed by 'box'
while ~box(prop_rnd)
  prop_rnd = mvnrnd(mu, sigma, N);
end

end
%% END