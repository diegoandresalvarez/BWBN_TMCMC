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

%%  External excitation
type_u = 1;         % Data Daniel
% type_u = 2;         % Data Kunnath
% type_u = 3;         % Data Gill
% type_u = 4;         % Data Chiu
% type_u = 5;         % Data Sharma
switch type_u
  case 1        % Data Daniel
    %% Define experimental data:
    m  = 456;                         % weight (kgf) = mass (kgm)
    % Remember that 1 kgm weighs 1 kgf.
    k  = 6.2684;                      % stiffness (kN/mm)

    % Here we multiply by 1000 because of the units. The data at the
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN_daniel.txt');
    u       = 1000*load_kN/m;
    N       = length(u);      % number of observations
    tt      = 0:0.02:((length(load_kN)-1)*0.02);
    tt      = tt';
    %% Real displacements (Measured)
    displ   = dlmread('data/displ_daniel.txt');
  case 2        % Data Kunnath
    %% Define experimental data:
    m  = 456;                         % weight (kgf) = mass (kgm)
    % Remember that 1 kgm weighs 1 kgf.
    k  = 6;                           % stiffness (kN/mm)

    % Here we multiply by 1000 because of the units. The data at the
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN_kunnath.txt');
    u       = 1000*load_kN/m;
    N       = length(u);      % number of observations
    tt      = 0:0.02:((length(load_kN)-1)*0.02);
    tt      = tt';
    %% Real displacements (Measured)
    displ   = dlmread('data/displ_kunnath.txt');
  case 3        % Data Gill
    %% Define experimental data:
    m  = 3200;                        % weight (kgf) = mass (kgm)
    % Remember that 1 kgm weighs 1 kgf.
    k  = 72.140;                      % stiffness (kN/mm)

    % Here we multiply by 1000 because of the units. The data at the
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN_gill.txt');
    u       = 1000*load_kN/m;
    N       = length(u);      % number of observations
    tt      = 0:0.02:((length(load_kN)-1)*0.02);
    tt      = tt';
    %% Real displacements (Measured)
    displ   = dlmread('data/displ_gill.txt');
  case 4        % Data Chiu
    %% Define experimental data:
    m  = 1360;                        % weight (kgf) = mass (kgm)
    % Remember that 1 kgm weighs 1 kgf.
    k  = 51.663;                      % stiffness (kN/mm)

    % Here we multiply by 1000 because of the units. The data at the
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN_chiu.txt');
    u       = 1000*load_kN/m;
    N       = length(u);      % number of observations
    tt      = 0:0.02:((length(load_kN)-1)*0.02);
    tt      = tt';
    %% Real displacements (Measured)
    displ   = dlmread('data/displ_chiu.txt');
  case 5        % Data Sharma
    %% Define experimental data:
    m  = 187.20;                      % weight (kgf) = mass (kgm)
    % Remember that 1 kgm weighs 1 kgf.
    k  = 5.8333;                      % stiffness (kN/mm)

    % Here we multiply by 1000 because of the units. The data at the
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN_sharma.txt');
    u       = 1000*load_kN/m;
    N       = length(u);      % number of observations
    tt      = 0:0.02:((length(load_kN)-1)*0.02);
    tt      = tt';
    %% Real displacements (Measured)
    displ   = dlmread('data/displ_sharma.txt');
  otherwise
    error('Invalid external excitation');
end

w0 = sqrt(k*(10^6)/m);              % natural frequency (rad/s)
dt = tt(2)-tt(1);                   % Runge-Kutta time step (sec)

%% Total dissipated energy (Measured form hysteresis cycles)
% Rememeber that the total dissipated energy is the area enclosed by the
% hysteresis cycles. The energy is a cummulative measure. The eq. is
% divided by 'm' (mass) because of the units, without this, the eq. has
% units of 'J'; with the factor, the eq. has units of 'J/kg'.
tot_diss_energy = cumtrapz(displ,load_kN)./m;
data            = [displ tot_diss_energy];

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
            1.85    % nu0             = strength degradation
            0.60    % deltanu         = strength degradation parameter
            2.20    % A0              = hysteresis amplitude
            0.78    % deltaA          = control parameter of 'A' with respect to the energy
            2.10    % eta0            = stiffness degradation
            1.15    % deltaeta        = stiffness degradation parameter

            % If you do not want to introduce Pinching effect, comment the
            % following six (6) lines. If you want to introduce Pinching
            % you must introduce also degradation
            5.00    % p               = parameter which controls initial pinching
            2.25    % vs0             = pinching severity
           -3.45    % psi0            = pinching parameter
            1.35    % deltapsi        = parameter which controls change of pinching
           -2.00    % lambda          = pinching parameter
            3.00    % q               = pinching parameter
        ]';

%% Setting boundaries of the parameters
%          lower    upper
%          bound    bound     comments
boundaries = [...
            eps      0.1          % sigma2_q_displ    Variance noise displacement
            eps      0.1          % sigma2_q_ene      Variance noise energy
            0.01     0.1          %  2  xi            damping ratio (Dimensionless)
            0.1      0.5          %  3  alpha         Bouc-Wen parameter
            0.5      4            %  4  beta          Bouc-Wen parameter
           -4        4            %  5  gamma         Bouc-Wen parameter
            1        3            %  6  n             Bouc-Wen parameter

            0.1      3            %  7  nu0           strength degradation
           -2        4            %  8  deltanu
            0.5      3            %  9  A0            strength and stiffness degradation
           -2        3            %  10 deltaA
            0.5      4            %  11 eta0          stiffness degradation
           -3        4            %  12 deltaeta
           
            0       10            %  13 p             pinching parameter
           -8        8            %  14 vs0           pinching parameter
           -4        4            %  15 psi0          pinching parameter
           -4        4            %  16 deltapsi      pinching parameter
           -4        4            %  17 lambda        pinching parameter
           -4        4            %  18 q             pinching parameter
         ]';

% Prior PDF
prior = @(x) prod(unifpdf(x, boundaries(1,:), boundaries(2,:)));

% Function that samples from the prior PDF
prior_rnd = @(NN) sample_unifrnd(boundaries, NN);

%% The log-likelihood of X given 'theta' - Likelihood PDF
log_f_X_theta = @(theta) log_p_X_theta(theta, data, u, w0, tt, dt, env);

%% TMCMC
%  Enable parallel computing within MATLAB
if matlabpool('size') == 0        % Activate all available cores
   matlabpool open
end;

% Start TMCMC
Nj = 100;                        % Number of samples per iteration
tic
[theta_fT_D, log_fD, p, Theta] = tmcmc(log_f_X_theta, prior, prior_rnd, Nj);
toc

matlabpool close                  % Deactivate all cores

%% Save results to a file
nc = length(Theta);

save('Results_BWBN.mat','Theta','p','boundaries','log_fD');
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