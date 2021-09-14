function plots_bwbn(x, t, param)
%% plots_bwbn(x, t, param)
%
% Plots Bouc-Wen-Baber-Noori
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            05-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------

%% Beginning
% u = normalized mass external excitation

npar = numel(param);

%           x(:,1)              % system displacement
xd        = x(:,2);             % system velocity
z         = x(:,3);             % hysteretic component
e         = x(:,4);             % dissipated energy
w0        = param(1);           % natural frequency
xi        = param(2);           % damping ratio
alpha     = param(3);           % alpha
beta      = param(4);           % beta
gamma     = param(5);           % gamma
n         = param(6);           % n

if  npar == 6
   nueps  = 1;                  % strength degradation
   Aeps   = 1;
   etaeps = 1;                  % stiffness degradation
   
   h = 1;                       % pinching function
end

if npar >= 7
   nu0      = param(7);      deltanu  = param(8);
   A0       = param(9);      deltaA   = param(10);
   eta0     = param(11);     deltaeta = param(12);
   
   nueps  = nu0  + deltanu *e;  % strength degradation
   Aeps   = A0   - deltaA  *e;
   etaeps = eta0 + deltaeta*e;  % stiffness degradation
end

if npar == 12
   h = 1;
end

if npar >= 13
   p        = param(13);
   vs0      = param(14);
   psi0     = param(15);
   deltapsi = param(16);
   lambda   = param(17);
   q        = param(18);
   
   zu  = (1./(nueps*(beta+gamma))).^(1/n);
   vs1 = (1 - exp(-p*e)).*vs0;
   vs2 = (psi0 + deltapsi.*e).*(lambda + vs1);
   h   = 1 - vs1.*exp(-((z.*sign(xd) - q.*zu).^2)./(vs2.^2));   % pinching degradation
end

%% Plots
% Plot C1 = h*A/eta
C1 = Aeps.*h./etaeps;
figure;
plot(t, C1);
title('A_{eps}*h/eta')

% Plot C2 = h*nu/eta
C2 = nueps;
figure;
plot(t, C2);
title('\nu_{eps}')

% Plot C2 = h*nu/eta
C3 = etaeps;
figure;
plot(t, C3);
title('\eta_{eps}')

% Plot C2 = h*nu/eta
C4 = h;
figure;
plot(t, C4);
title('h')

% Plot C2 = h*nu/eta
figure;
plot(t, zu);
title('z_{ueps}')

end
%% END