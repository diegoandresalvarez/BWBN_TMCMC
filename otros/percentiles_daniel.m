%% compute percentiles of data

%% Load real data
real_displ  = load('Results/Results_Daniel_PC_abajo_26_Oct_2013/plots/displ_real.txt');
real_energy = load('Results/Results_Daniel_PC_abajo_26_Oct_2013/plots/energy_real.txt');
t           = load('Results/Results_Daniel_PC_abajo_26_Oct_2013/plots/time.txt');
% t           = load('/home/gilberto/Documents/Bayesian_estimation/programs/tmcmc/Paper_computers_and_structures/fig/fourier/time.txt');

%% Load parameters
load('Results/Results_Daniel_PC_abajo_26_Oct_2013/Results_BWBN.mat');
params = Theta{end}(:,3:end)';

%% load estimated displacements (those with maximum probability)
displ_est  = load('Results/Results_Daniel_PC_abajo_26_Oct_2013/plots/displ_est.txt');
energy_est = load('Results/Results_Daniel_PC_abajo_26_Oct_2013/plots/energy_est.txt');
% displ_est  = load('/home/gilberto/Documents/Bayesian_estimation/programs/tmcmc/Paper_computers_and_structures/fig/fourier/displ_est_sim1.txt');
% energy_est = load('/home/gilberto/Documents/Bayesian_estimation/programs/tmcmc/Paper_computers_and_structures/fig/fourier/energy_est.txt');

%% Define experimental data:
m  = 456;                         % weight (kgf) = mass (kgm)
% Remember that 1 kgm weighs 1 kgf.
k  = 6.2684;                      % stiffness (kN/mm)
% Experimental load
type_u = 1;
% type_u = 2;

l = length(real_displ);
% l = length(displ_est);
L = size(params,2);

displ           = zeros(l,L);
Fz              = zeros(l,L);
tot_diss_energy = zeros(l,L);

%% Compute BWBN
for i = 1:L
  [displ(:,i), Fz(:,i), tot_diss_energy(:,i)] = bwbn_function(m,k,type_u,params(:,i));
end

%% compute mean displacement and mean energy
mean_displ  = mean(displ,2);
mean_energy = mean(tot_diss_energy,2);

%% Calculate the 25th, 50th, and 75th percentiles of the displacements and
%  the energy
per_displ  = prctile(displ,[10 25 50 75 90],2);
per_energy = prctile(tot_diss_energy,[10 25 50 75 90],2);

%% Compute interdecile range
range_displ  = per_displ(:,5)  - per_displ(:,1);
range_energy = per_energy(:,5) - per_energy(:,1);

%{
%% save results
dlmwrite('percentil_10_25_50_75_90_displ.txt',per_displ);
dlmwrite('percentil_10_25_50_75_90_energy.txt',per_energy);
%}

%% Plots
% 1. Displacements
figure
subplot(2,1,1)
  fill([t; flipud(t)],[per_displ(:,1); flipud(per_displ(:,5))],[14 14 14]/20, 'EdgeColor', 'r');
  hold on
  fill([t; flipud(t)],[per_displ(:,2); flipud(per_displ(:,4))],[14 14 14]/15, 'EdgeColor', 'm');
  plot(t,real_displ,'b');
  plot(t,displ_est,'k','lineWidth',2);
  plot(t,mean_displ,'--g','lineWidth',2);
  plot(t,per_displ(:,3),'c');
  xlabel('Time (s)', 'FontSize', 16);
  ylabel('Displacement (mm)', 'FontSize', 16);
  title('Time vs. Displacement', 'FontSize', 18);
  legend('p_{10}-p_{90}','p_{25}-p_{75}','Real','MAP Estimated','Mean','p_{50}','location','Best');
%   legend('p_{10}-p_{90}','p_{25}-p_{75}','MAP Estimated','Mean','p_{50}','location','Best');
  grid on
subplot(2,1,2)
  plot(t,range_displ);
  xlabel('Time (s)', 'FontSize', 16);
  ylabel('Interdecile range (mm)', 'FontSize', 16);
  grid on

% 2. Total dissipated energy
figure
subplot(2,1,1)
  fill([t; flipud(t)],[per_energy(:,1); flipud(per_energy(:,5))],[14 14 14]/20, 'EdgeColor', 'r');
  hold on
  fill([t; flipud(t)],[per_energy(:,2); flipud(per_energy(:,4))],[14 14 14]/15, 'EdgeColor', 'm');
  plot(t,real_energy,'b');
  plot(t,energy_est,'k','lineWidth',2);
  plot(t,mean_energy,'--g','lineWidth',2);
  plot(t,per_energy(:,3),'c');
  xlabel('Time (s)', 'FontSize', 16);
  ylabel('Total dissipated energy (J/kg)', 'FontSize', 16);
  title('Time vs. Total dissipated energy', 'FontSize', 18);
  legend('p_{10}-p_{90}','p_{25}-p_{75}','Real','MAP Estimated','Mean','p_{50}','location','Best');
%   legend('p_{10}-p_{90}','p_{25}-p_{75}','MAP Estimated','Mean','p_{50}','location','Best');
  grid on
subplot(2,1,2)
  plot(t,range_energy);
  xlabel('Time (s)', 'FontSize', 16);
  ylabel('Interdecile range (J/kg)', 'FontSize', 16);
  grid on