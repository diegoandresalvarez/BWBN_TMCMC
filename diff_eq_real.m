function dxdt = diff_eq_real(x,u,param)
%%  dxdt = diff_eq_real(x,u,param)
%   This function returns the set of equations of the Bouc-Wen hysteresis
%   model with degradation and pinching.
%
%   Bibliography:
%
%   - FOLIENTE, Greg C. "Hysteresis modeling of wood joints and structural
%     systems". Journal of Structural Engineering. Vol. 121. Nro. 6. June.
%     1995.
%
% -------------------------------------------------------
% | Developed by:   Diego Andres Alvarez Marin          |
% |                 diegotorquemada@gmail.com           |
% |                 National University of Colombia     |
% |                 Manizales, Colombia.                |
% |                                                     |
% |                 Gilberto Alejandro Ortiz Garcia     |
% |                 gialorga@gmail.com                  |
% |                 National University of Colombia     |
% |                 Manizales, Colombia.                |
% -------------------------------------------------------
%
%   Date: 25 - Aug - 2011

%% Beginning:
% u = normalized mass external excitation

npar = numel(param);

%           x(1)                % system displacement
xd        = x(2);               % system velocity
z         = x(3);               % hysteretic component
e         = x(4);               % dissipated energy
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
   
   zu  = (1/(nueps*(beta+gamma)))^(1/n);
   vs1 = (1 - exp(-p*e))*vs0;
   vs2 = (psi0 + deltapsi*e)*(lambda + vs1);
   h   = 1 - vs1*exp(-((z*sign(xd) - q*zu)^2)/(vs2^2));   % pinching degradation
end

% NOTES:
% - In the second eq., 'u' is multiplied by 1000 because of the units, so
%   'mm/s^2' will be obtained.
%
% - The last eq. was multiplied by 1e-6 because of the units, without the
%   factor, the eq. has units of 'mm^2/s^2'; with the factor, the last eq.
%   (dissipated hysteretic energy) has units of 'J/kg'.
dxdt = [...
        xd;
        1000*u - 2*xi*w0*xd - alpha*w0^2*x(1) - (1-alpha)*w0^2*z;
        h*(Aeps*xd - nueps*(beta*abs(xd)*(abs(z)^(n-1))*z + gamma*xd*(abs(z)^n)))/etaeps;
        1e-6*(1-alpha)*w0^2*z*xd
       ];

% In case the Bouc-Model of hysteresis diverges, we set the response to Inf
if any(any(isnan(dxdt))) || any(any(isinf(dxdt)))
 dxdt = inf(4,1);
end;

end
%% END