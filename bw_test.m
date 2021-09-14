%% Bouc-Wen-Baber-Noori (BWBN) hysteresis model
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
% -------------------------------------------------------
% | Developed by:   Diego Andres Alvarez Marin          |
% |                 diegotorquemada@gmail.com           |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% |                                                     |
% |                 Gilberto Alejandro Ortiz Garcia     |
% |                 gialorga@gmail.com                  |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% -------------------------------------------------------
%
%   Date: 09 - Sep - 2011

%% Beginning:
clear, clc, close all
tic
%% loading experimental data:
m  = 456;                           % weight (kgf) = mass (kgm)
% Remember that 1 kgm weighs 1 kgf.
k  = 6.2684;                        % stiffness (kN/mm)
w0 = sqrt(k*(10^6)/m);              % natural frequency (rad/s)

%%  External excitation
type_u = 1;                         % Force (kN)
% type_u = 2;                         % Acceleration (mm/s^2)
% type_u = 3;                         % sinusoidal;
% type_u = 4;                         % Pattern ASTM-E2126 (Load (kN))
% type_u = 5;                         % White noise (Random excitation)
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
%     u_m = dlmread('data/japan_strong.txt');
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

dt      = tt(2)-tt(1);                    % Runge-Kutta time step (sec)

%% Setting Bouc-Wen parameters:
xi        =   0.0878;
alpha     =   0.1997;
beta      =   1.9420;
gamma     =  -0.4643;
n         =   1.0512;

% If you do not want to introduce degradation effect, comment the following
% six (6) lines, and also comment the section including 'pinching' in such
% a case.
nu0       =   0.1397;
deltanu   =   0.9759;
A0        =   0.5001;
deltaA    =   0.0355;
eta0      =   0.6750;
deltaeta  =   3.9255;

% If you do not want to introduce Pinching effect, comment the following
% six (6) lines. If you want to introduce Pinching you must introduce also
% degradation
p         =   3.6874;
vs0       =   6.2891;
psi0      =   1.5485;
deltapsi  =  -3.3111;
lambda    =  -3.6275;
q         =  -1.3256;

param = [...
            w0            % natural frequency (rad/s)
            xi           % damping ratio (Adim.)
            alpha        % ratio of post-yield to pre-yield stiffness
            beta         % Bouc-Wen parameter
            gamma        % Bouc-Wen parameter
            n            % Bouc-Wen parameter

            % If you do not want to introduce degradation effect, comment the
            % following six (6) lines, and also comment the section including 
            % 'pinching' in such a case.
            nu0          % strength degradation
            deltanu      % strength degradation parameter
            A0           % hysteresis amplitude
            deltaA       % control parameter of 'A' with respect to the energy
            eta0         % stiffness degradation
            deltaeta     % stiffness degradation parameter

            % If you do not want to introduce Pinching effect, comment the
            % following six (6) lines. If you want to introduce Pinching you
            % must introduce also degradation
            p            % parameter which controls initial pinching
            vs0          % pinching severity
            psi0         % pinching parameter
            deltapsi     % parameter which controls change of pinching
            lambda       % pinching parameter
            q            % pinching parameter
        ];

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
x_k = x_k';

%% Computing "Restoring Force" (Fz):
Fz = alpha*k*x_k(:,1) + (1-alpha)*k*x_k(:,3);
toc

%% Plot the results:

%% Displacement:
figure;
plot(tt,x_k(:,1),'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Displacement (mm)', 'FontSize', 16);
title('Time vs. Displacement', 'FontSize', 18);
grid on

%% Hysteresis cycles:
figure;
plot(x_k(:,1),Fz,'b');
xlabel('Displacement (mm)', 'FontSize', 16);
ylabel('Restoring force (kN)', 'FontSize', 16);
title('Displacement vs. Restoring force', 'FontSize', 18);
grid on

%% Total dissipated energy
% Compute 'dissipated elastic energy'
%
%                                            /t_f
%  dissipated_elastic_energy = alpha * w0^2 *|     est_displ * est_vel dt
%                                            /t_0
%
% Rememeber that this is a cummulative measure.
%
diss_elastic_energy = 1e-6*alpha*(w0^2)*cumtrapz(tt, x_k(:,1).*x_k(:,2));

% tot_diss_energy = diss_elastic_energy + diss_hysteretic_energy
tot_diss_energy = diss_elastic_energy + x_k(:,4);

figure;
plot(tt,tot_diss_energy,'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Dissipative energy ( J / kg )', 'FontSize', 16);
title('Time vs. Dissipative energy', 'FontSize', 18);
grid on

%% END