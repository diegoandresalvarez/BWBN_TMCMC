function prior_rnd = sample_unifrnd(boundaries, N)
%% prior_rnd = sample_unifrnd(boundaries, N)
%
% Get 'N' random samples from a continuous uniform distribution on the
% interval from boundaries(1,:) to boundaries(2,:).
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            11-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning
n         = size(boundaries,2);
prior_rnd = zeros(N,n);

for i = 1:N
  prior_rnd(i,:) = unifrnd(boundaries(1,:), boundaries(2,:));
end
end
%% END