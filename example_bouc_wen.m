%% EXAMPLE
%
%               x_k = F(x_{k-1}, u_{k-1}) + v_{k-1}     (1)
%               z_k = H(x_k, u_k)         + w_k         (2)
%
% BIBLIOGRAPHY:
%
% - FOLIENTE, Greg C. "Hysteresis modeling of wood joints and structural
%   systems". Journal of Structural Engineering. Vol. 121. Nro. 6. June.
%   1995.
%
% - CHING, Jianye; CHEN, Yi-Chun. "Transitional Markov Chain Monte Carlo
%   method for Bayesian model updating, model class selection and model
%   averaging". Journal of Engineering Mechanics. ASCE. 133(7):816-832,
%   July 1, 2007.
%
% - MUTO, Matthew; BECK, James L. "Bayesian updating and model class
%   selection for hysteretic structural models using stochastic
%   simulation". Journal of vibration and control. 14(1-2):7-34, 2008.
%
% -------------------------------------------------------------------------
% * Developed by:                Date:            Mail:
%   Gilberto A. Ortiz            11-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning
close all; clear all; clc;

%% Store seed for random number generator
seed = rng;

%% Define experimental data:
m  = 500;                           % weight (kgf) = mass (kgm)
% Remember that 1 kgm weighs 1 kgf.
k  = 6;                             % stiffness (kN/mm)
w0 = sqrt(k*(10^6)/m);              % natural frequency (rad/s)

%%  External excitation
% type_u = 1                        % Force (kN)
% type_u = 2                        % Acceleration (mm/s^2)
% type_u = 3                        % Acceleration (mm/s^2)
type_u = 4;                         % sinusoidal
% type_u = 5                        % Pattern ASTM-E2126 (Load (kN))
% type_u = 6                        % White noise (Random excitation)
switch type_u
  case 1        % Force (kN)
    % Here we multiply by 1000 because of the units. The data at the
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN.txt');
    u       = 1000*load_kN/m;
    N       = length(u);      % number of observations
    tt      = 0:0.02:(length(load_kN)*0.02);
  case 2        % Acceleration (mm/s^2)
    % Here we multiply by 9.81 because of the units. The earthquake has
    % units of '%g', and we need acceleration in 'm/s^2'.
    u_m = dlmread('data/Northridge_1994_25_sec.txt');
    tt  = u_m(:,1)';
    u   = 9.81*u_m(:,2)';
    N   = length(u);          % number of observations
  case 3        % Acceleration (mm/s^2);
    % Here we multiply by 9.81 because of the units. The earthquake has
    % units of '%g', and we need acceleration in 'm/s^2'.
    u_m = dlmread('data/japan_2011.txt');
    tt  = u_m(:,1)';
    u   = 9.81*u_m(:,2)';
    N   = length(u);          % number of observations
  case 4        % Sinusoidal
    % Simulated input signal (sinusoidal case)
    tt  = (0:0.02:19.98)';
    uu  = tt.*sin(2*pi*tt);   % kN
    u   = 1000*uu/m;          % m/s^2
    N   = length(u);          % number of observations
  case 5        % Pattern ASTM-E2126 (Load (kN))
    % Here we multiply by 1000 because of the units. The data returned by
    % function 'generate_signal_astm' is in 'kN', and we need acceleration
    % in 'm/s^2'.
    [tt, u_m] = generate_signal_astm(20, 0.02);
    u         = 1000*u_m'/m;
    N         = length(u);    % number of observations
  case 6        % White noise (Random excitation)
    N  = 1250;                % number of observations
    tt = linspace(0,25,N);    % time vector
    u  = randn(1,N);
  otherwise
    error('Invalid external excitation');
end

dt = tt(2)-tt(1);                 % Runge-Kutta time step (sec)

%% Variance addititve white noise (displacements)
sigma2_q_disp = 0.01;

%% Setting Bouc-Wen parameters:
xi        =   0.05;
alpha     =   0.20;
beta      =   2.00;
gamma     =  -1.00;
n         =   1.20;

% If you do not want to introduce degradation effect, comment the following
% six (6) lines, and also comment the section including 'pinching' in such
% a case.
% nu0       =   1.20;
% deltanu   =   0.40;
% A0        =   1.10;
% deltaA    =   0.10;
% eta0      =   1.30;
% deltaeta  =   0.50;

% If you do not want to introduce Pinching effect, comment the following
% six (6) lines. If you want to introduce Pinching you must introduce also
% degradation
% p         =   2.00;
% vs0       =   1.00;
% psi0      =   0.50;
% deltapsi  =   0.60;
% lambda    =   0.90;
% q         =   1.00;

param = [...
            w0        % natural frequency (rad/s)
            xi        % damping ratio (Adim.)
            alpha     % ratio of post-yield to pre-yield stiffness
            beta      % Bouc-Wen parameter
            gamma     % Bouc-Wen parameter
            n         % Bouc-Wen parameter

            % If you do not want to introduce degradation effect, comment
            % the following six (6) lines, and also comment the section
            % including 'pinching' in such a case.
%            nu0       % strength degradation
%            deltanu   % strength degradation parameter
%            A0        % hysteresis amplitude
%            deltaA    % control parameter of 'A' with respect to the energy
%            eta0      % stiffness degradation
%            deltaeta  % stiffness degradation parameter

            % If you do not want to introduce Pinching effect, comment the
            % following six (6) lines. If you want to introduce Pinching
            % you must introduce also degradation
%            p         % parameter which controls initial pinching
%            vs0       % pinching severity
%            psi0      % pinching parameter
%            deltapsi  % parameter which controls change of pinching
%            lambda    % pinching parameter
%            q         % pinching parameter
        ]';

%% Definition of the state function
BW_real = @(x,u) diff_eq_real(x,u,param);
Fexact  = @(x_k,u_k) rk_discrete(BW_real,x_k,u_k,dt);

%% Initial condition:
% Initial condition x, xd, z, e (Displacement, velocity. hysteretic
% displacement, hysteretic energy) set to zero
x_0 = zeros(4,1);

% vector where I'm going to store the system response
x_k = zeros(4,N);

x_k(:,1) = x_0;           % Initial state

%% Computing system response
for i = 2:N
  x_k(:,i) = Fexact(x_k(:,i-1),u(i-1));
end
x_k = x_k';

displ_noisy = x_k(:,1) + sqrt(sigma2_q_disp)*randn(length(x_k(:,1)),1);
displ_noisy(1,1) = 0;

%% Computing "Restoring Force" (Fz):
Fz = alpha*k*x_k(:,1) + (1-alpha)*k*x_k(:,3);

%% Total dissipated energy
% Compute 'dissipated elastic energy'
%
%                                            /t_f
%  dissipated_elastic_energy = alpha * w0^2 *|     est_displ * est_vel dt
%                                            /t_0
%
% Rememeber that this is a cummulative measure. The eq. is multiplied by
% 1e-6 because of the units, without the factor, the eq. has units of
% 'mm^2/s^2'; with the factor, the eq. (dissipated hysteretic energy) has
% units of 'J/kg'.

diss_elastic_energy = 1e-6*alpha*(w0^2)*cumtrapz(tt, displ_noisy.*x_k(:,2));

% tot_diss_energy = diss_elastic_energy + diss_hysteretic_energy
tot_diss_energy = diss_elastic_energy + x_k(:,4);

%{
%% Plot the results:

%% Displacement:
figure;
subplot(2,2,1);
plot(tt,x_k(:,1),'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Displacement (mm)', 'FontSize', 16);
title('Time vs. Displacement (Simulated system)', 'FontSize', 16);
grid on

subplot(2,2,3);
plot(tt,tot_diss_energy,'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Dissipative energy ( J / kg )', 'FontSize', 16);
title('Time vs. Dissipative energy (Simulated system)', 'FontSize', 16);
grid on

%% Hysteresis cycles:
subplot(2,2,[2 4]);
plot(x_k(:,1),Fz,'b');
xlabel('Displacement (mm)', 'FontSize', 16);
ylabel('Restoring force (kN)', 'FontSize', 16);
title('Displacement vs. Restoring force (Simulated system)', 'FontSize', 16);
grid on
%}

%% Real displacements and dissipated energy (Simulated)
data = [displ_noisy tot_diss_energy];

%% Compute envelopes for normalization of data between [-1, 1]

% Envelope displacements
der  = diff(data);
idx1 = find(der(1:end-1,1).*der(2:end,1)<=0)+1;
idx2 = find(der(1:end-1,2).*der(2:end,2)<=0)+1;
env1 = interp1([0; tt(idx1)],[0; abs(data(idx1,1))],tt,'linear','extrap');
env2 = interp1([0; tt(idx2)],[0; abs(data(idx2,2))],tt,'linear','extrap');

env = [env1 env2];

% change zero for 'eps'
env(env == 0) = eps;

%{
figure
hold on
plot(tt,data(:,1),'bx-', tt(idx1),data(idx1,1),'ro', tt,env(:,1),'g');

figure
hold on
plot(tt,data(:,2),'bx-', tt(idx2),data(idx2,2),'ro', tt,env(:,2),'g');

figure
plot(tt,data(:,1)./env(:,1));

figure
plot(tt,data(:,2)./env(:,2));

dlmwrite('weighting.txt', env, 'precision', '%.6f', 'roffset', 1, 'delimiter', '\n');
%}

%% ------------ Transitional Markov Chain Monte Carlo (TMCMC) ------------

%% Define Prior hyperparameters and prior PDF
theta_0 = [...
            0.10    % sigma2_q_displ  = Variance noise displacement
            0.01    % sigma2_q_ene    = Variance noise energy
            0.15    % xi              = damping ratio (Adim.)
            0.15    % alpha           = ratio of post-yield to pre-yield stiffness
            1.50    % beta            = Bouc-Wen parameter
           -2.00    % gamma           = Bouc-Wen parameter
            1.00    % n               = Bouc-Wen parameter

            % If you do not want to introduce degradation effect, comment
            % the following six (6) lines, and also comment the section
            % including 'pinching' in such a case.
%             1.85    % nu0              = strength degradation
%             0.60    % deltanu          = strength degradation parameter
%             2.20    % A0               = hysteresis amplitude
%             0.78    % deltaA           = control parameter of 'A' with respect to the energy
%             2.10    % eta0             = stiffness degradation
%             1.15    % deltaeta         = stiffness degradation parameter

            % If you do not want to introduce Pinching effect, comment the
            % following six (6) lines. If you want to introduce Pinching
            % you must introduce also degradation
%             5.00    % p                = parameter which controls initial pinching
%             2.25    % vs0              = pinching severity
%            -3.45    % psi0             = pinching parameter
%             1.35    % deltapsi         = parameter which controls change of pinching
%            -2.00    % lambda           = pinching parameter
%             3.00    % q                = pinching parameter
        ]';

%% Setting boundaries of the parameters
%          lower    upper
%          bound    bound     comments
boundaries = [...
            eps      0.1          % sigma2_q_displ    Variance noise displacement
            eps      0.1          % sigma2_q_ene      Variance noise energy
            0.01     0.1          %  2  xi            damping ratio (Dimensionless)
            0.1      0.3          %  3  alpha         Bouc-Wen parameter
            1        4            %  4  beta          Bouc-Wen parameter
           -4        4            %  5  gamma         Bouc-Wen parameter
            1        3            %  6  n             Bouc-Wen parameter

%             0.1      3            %  7  nu0           strength degradation
%            -2        4            %  8  deltanu
%             0.5      3            %  9  A0            strength and stiffness degradation
%            -2        3            %  10 deltaA
%             0.5      4            %  11 eta0          stiffness degradation
%            -3        4            %  12 deltaeta

%             0       10            %  13 p             pinching parameter
%            -8        8            %  14 vs0           pinching parameter
%            -4        4            %  15 psi0          pinching parameter
%            -4        4            %  16 deltapsi      pinching parameter
%            -4        4            %  17 lambda        pinching parameter
%            -4        4            %  18 q             pinching parameter
         ]';

% Prior PDF
prior = @(x) prod(unifpdf(x, boundaries(1,:), boundaries(2,:)));

% Function that samples from the prior PDF
prior_rnd = @(NN) sample_unifrnd(boundaries, NN);

%% The log-likelihood of X given 'theta' - Likelihood PDF
log_f_X_theta = @(theta) log_p_X_theta(theta, data, u, w0, tt, dt, env);

%%  Activate all available cores (parallel computing) within MATLAB
% parpool;

%% Start TMCMC
Nj = 1000;                       % Number of samples per iteration
burnin = 50;                     % parameter of the Metropolis-Hastings
last_burnin = 200;               % burn-in in the last iteration

tic
[theta_fT_D, log_fD, p, Theta] = tmcmc(log_f_X_theta, prior, prior_rnd, Nj, burnin, last_burnin);
toc

%% Deactivate all cores
% delete(gcp('nocreate'));          % Deactivate all cores

%% Save results to a file
nc = length(Theta);

save('Results_BW.mat','Theta','p','boundaries','log_fD');
%{
for mm = 1:nc
  dlmwrite('Results_BW.txt',Theta{mm},'-append');
  dlmwrite('Results_BW.txt',' ','delimiter','\n','-append');
end
%}

%{
%% Plot TMCMC results
np     = size(Theta{end},2);
Labels = {'\sigma_{q_{dis}}','\sigma_{q_{ene}}','\xi','\alpha','\beta',...
          '\gamma','n','\nu_{0}','\delta_{\nu}','A_{0}','\delta_{A}',...
          '\eta_{0}','\delta_{\eta}','p','\zeta_{0}','\psi_{0}',...
          '\delta_{\psi}','\lambda','q'};

for j = 1:nc
  figure(j);
  [H,AX,BigAx,P]=plotmatrix(Theta{j},'.b');
  title(sprintf('TMCMC results - Bouc Wen model - j = %d', j), 'FontSize', 18);

  for i = 1:np
    set(get(AX(i,1),'ylabel'),'string',Labels(i),'FontSize', 15);
    set(get(AX(np,i),'xlabel'),'string',Labels(i),'FontSize', 15);
    set(P(i),'EdgeColor','k','FaceColor','c');
  end
end
%}

%% END