function [displ, Fz, tot_diss_energy] = bwbn_function(m,k,type_u,params)
%% Bouc-Wen-Baber-Noori (BWBN) hysteresis model
%
% [displ, Fz, tot_diss_energy] = bwbn_function(m,k,type_u,params)
%
% xdd + 2*xi*w0*xd + alpha*w0^2*x(1) + (1-alpha)*w0^2*z = u;
% 
% nueps  = nu0  + deltanu *e;
% Aeps   = A0   - deltaA  *e;
% etaeps = eta0 + deltaeta*e;
%
% zu  = (1/(nueps*(beta+gamma)))^(1/n);
% vs1 = (1 - exp(-p*e))*vs0;
% vs2 = (psi0 + deltapsi*e)*(lambda + vs1);
% h   = 1 - vs1*exp(-((z*sign(xd) - q*zu)^2)/(vs2^2));
%
%
%      h*(Aeps*xd - nueps*(beta*abs(xd)*(abs(z)^(n-1))*z + gamma*xd*(abs(z)^n)))
% zd = -------------------------------------------------------------------------
%                               etaeps;
%
%   (*) Input data:
%   - m:      weight (kgf) = mass (kgm)
%   - k:      stiffness (kN/mm)
%   - type_u: 1 = Force (kN)
%             2 = Acceleration (mm/s^2)
%             3 = Sinusoidal;
%             4 = Pattern ASTM-E2126 (Load (kN))
%             5 = White noise (Random excitation)
%   - params: BWBN parameters:
%
%      params = [...
%            xi           % damping ratio (Adim.)
%            alpha        % ratio of post-yield to pre-yield stiffness
%            beta         % Bouc-Wen parameter
%            gamma        % Bouc-Wen parameter
%            n            % Bouc-Wen parameter
%
%            % If you do not want to introduce degradation effect, comment
%            % the following six (6) lines, and also comment the section
%            % including 'pinching' in such a case.
%            nu0          % strength degradation
%            deltanu      % strength degradation parameter
%            A0           % hysteresis amplitude
%            deltaA       % control parameter of 'A' with respect to the energy
%            eta0         % stiffness degradation
%            deltaeta     % stiffness degradation parameter

%            % If you do not want to introduce Pinching effect, comment the
%            % following six (6) lines. If you want to introduce Pinching
%            % you must introduce also degradation
%            p            % parameter which controls initial pinching
%            vs0          % pinching severity
%            psi0         % pinching parameter
%            deltapsi     % parameter which controls change of pinching
%            lambda       % pinching parameter
%            q            % pinching parameter
%        ];
%
%   (*) Output data:
%   - displ:           Displacement (mm)
%   - Fz:              Restoring force (kN)
%   - tot_diss_energy: Total dissipated energy (J/kg)
%
%   Bibliography:
%
%   - FOLIENTE, Greg C. "Hysteresis modeling of wood joints and structural
%     systems". Journal of Structural Engineering. Vol. 121. Nro. 6. June.
%     1995.
%
%   - ASTM E2126-11: "Standard test methods for cyclic (reversed) load test for
%     shear resistance of vertical elements of the lateral force resisting
%     systems for buildings".
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            05-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------

%% Beginning:

%% natural frequency (rad/s):
w0 = sqrt(k*(10^6)/m);              % 

%%  External excitation
switch type_u
  case 1        % Force (kN);   
    % Here we multiply by 1000 because of the units. The data at the
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN_daniel.txt');
    u       = 1000*load_kN/m;
    N       = length(u);      % number of observations
    tt      = 0:0.02:((length(load_kN)-1)*0.02);
  case 2        % Acceleration (mm/s^2);
    % Here we multiply by 9.81 because of the units. The earthquake has
    % units of '%g', and we need acceleration in 'm/s^2'.
    u_m = dlmread('data/Northridge_1994_25_sec.txt');
    tt  = u_m(:,1)';
    u   = 9.81*u_m(:,2)';
    N   = length(u);          % number of observations
  case 3        % Sinusoidal
    % Simulated input signal (sinusoidal case)
    tt  = (0:0.02:19.98)';
    uu  = tt.*sin(2*pi*tt);   % kN
    u   = 1000*uu/m;          % m/s^2
    N   = length(u);          % number of observations
  case 4        % Pattern ASTM-E2126 (Load (kN))
    % Here we multiply by 1000 because of the units. The data returned by
    % function 'generate_signal_astm' is in 'kN', and we need acceleration
    % in 'm/s^2'.
    [tt, u_m] = generate_signal_astm(20, 0.02);
    u         = 1000*u_m'/m;
    N         = length(u);    % number of observations
  case 5        % White noise (Random excitation)
    % White noise
    N  = 1250;                    % number of observations
    tt = linspace(0,25,N);      % time vector
    u  = randn(1,N);
  otherwise
    error('Invalid external excitation');
end

dt    = tt(2)-tt(1);                    % Runge-Kutta time step (sec)
param = [w0; params];

%% Definition of the state function
BW_real = @(x,u) diff_eq_real(x,u,param);
Fexact  = @(x_k,u_k) rk_discrete(BW_real,x_k,u_k,dt);

%% Initial condition:
% Initial condition x, xd, z, e (Displacement, velocity. hysteretic
% displacement, hysteretic energy) set to zero
x_0 = zeros(4,1);
x_k = zeros(length(x_0)); % vector where I'm going to store the system response

x_k(:,1) = x_0;           % Initial state

%% Computing system response
for i = 2:N
  x_k(:,i) = Fexact(x_k(:,i-1),u(i-1));
end
x_k   = x_k';
displ = x_k(:,1);

%% Computing "Restoring Force" (Fz):
alpha = params(2);
Fz    = alpha*k*displ + (1-alpha)*k*x_k(:,3);

%% Total dissipated energy
% Compute 'dissipated elastic energy'
%
%                                            /t_f
%  dissipated_elastic_energy = alpha * w0^2 *|     est_displ * est_vel dt
%                                            /t_0
%
% Rememeber that this is a cummulative measure.
%
diss_elastic_energy = 1e-6*alpha*(w0^2)*cumtrapz(tt, displ.*x_k(:,2));

% tot_diss_energy = diss_elastic_energy + diss_hysteretic_energy
tot_diss_energy = diss_elastic_energy + x_k(:,4);

%{
%% Plot the results:

%% Displacement:
figure;
plot(tt,displ,'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Displacement (mm)', 'FontSize', 16);
title('Time vs. Displacement', 'FontSize', 18);
grid on

%% Hysteresis cycles:
figure;
plot(displ,Fz,'b');
xlabel('Displacement (mm)', 'FontSize', 16);
ylabel('Restoring force (kN)', 'FontSize', 16);
title('Displacement vs. Restoring force', 'FontSize', 18);
grid on

figure;
plot(tt,tot_diss_energy,'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Dissipative energy ( J / kg )', 'FontSize', 16);
title('Time vs. Dissipative energy', 'FontSize', 18);
grid on
%}

end
%% END