function plot_TMCMC(Theta,boundaries,col1,col2,col3)
%% Routine to plot evolution of parameters computed with TMCM
%
%     (*) Input data:
%
%     - Theta:      Evolution of parameters (1 x nc cell array)
%     - boundaries: Boundaries of parameters
%     - col1, col2: Columns of Theta{i} to plot
%     - col3:       Third column of Theta to plot (optional)
%
%     (*) Output data
%
%     - Plots with the evolution of parameters
%
%  Modify the algorithm to suit your needs. Plot the parameters you want
%  and edit the labels accordingly.
%
% BIBLIOGRAPHY:
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
%   Gilberto A. Ortiz            04-Sep-2013      gialorga@gmail.com
%
%   Universidad Nacional de Colombia at Manizales. Civil Eng. Dept.
% -------------------------------------------------------------------------
%
%% Beginning

%% Compute size of cell array
[~, nc] = size(Theta);
N       = size(Theta{1},1);         % Number of elements

Labels = {'\sigma_{q_{dis}}','\sigma_{q_{ene}}','\xi','\alpha','\beta',...
          '\gamma','n','\nu_{0}','\delta_{\nu}','A_{0}','\delta_{A}',...
          '\eta_{0}','\delta_{\eta}','p','\zeta_{0}','\psi_{0}',...
          '\delta_{\psi}','\lambda','q'};

if nargin == 4
  for j = 1:nc
    figure(j)
      subplot(2,2,2);
        hist(Theta{j}(:,col1), ceil(sqrt(N)));
        xlabel(Labels(col1),'FontSize',16);
        ylabel('Frequency','FontSize',16);
        title('Histograms','FontSize',18);
      subplot(2,2,4);
        hist(Theta{j}(:,col2), ceil(sqrt(N)));
        xlabel(Labels(col2),'FontSize',16);
        ylabel('Frequency','FontSize',16);
      subplot(2,2,[1 3]);
        plot(Theta{j}(:,col1), Theta{j}(:,col2), 'r.');
        axis([boundaries(1,col1) boundaries(2,col1) boundaries(1,col2) boundaries(2,col2)])
        xlabel(Labels(col1),'FontSize',16);
        ylabel(Labels(col2),'FontSize',16);
      grid on;
      title(sprintf('Samples of f_{%d}', j-1),'FontSize',18);
  end

  figure
  hold on
  col = lines(nc);
  legendInfo = cell(1,nc);
  for j = 1:nc
    plot(Theta{j}(:,col1), Theta{j}(:,col2), '.', 'color', col(j,:));
    xlabel(Labels(col1),'FontSize',16);
    ylabel(Labels(col2),'FontSize',16);
    grid on;
    title('Evolution of parameters','FontSize',18);
    legendInfo{j} = sprintf('Stage = %d', j-1);
  end
  legend(legendInfo);
elseif nargin == 5
  for j = 1:nc
    figure(j)
      subplot(2,3,1);
        hist(Theta{j}(:,col1), ceil(sqrt(N)));
        xlabel(Labels(col1),'FontSize',16);
        ylabel('Frequency','FontSize',16);
      subplot(2,3,2);
        hist(Theta{j}(:,col2), ceil(sqrt(N)));
        xlabel(Labels(col2),'FontSize',16);
        ylabel('Frequency','FontSize',16);
        title('Histograms','FontSize',18);
      subplot(2,3,3);
        hist(Theta{j}(:,col3), ceil(sqrt(N)));
        xlabel(Labels(col3),'FontSize',16);
        ylabel('Frequency','FontSize',16);
      subplot(2,3,[4 6]);
        plot3(Theta{j}(:,col1), Theta{j}(:,col2), Theta{j}(:,col3), 'r.');
        axis([boundaries(1,col1) boundaries(2,col1)...
              boundaries(1,col2) boundaries(2,col2)...
              boundaries(1,col3) boundaries(2,col3)]);
        xlabel(Labels(col1),'FontSize',16);
        ylabel(Labels(col2),'FontSize',16);
        zlabel(Labels(col3),'FontSize',16);
      grid on;
      title(sprintf('Samples of f_{%d}', j-1),'FontSize',18);
  end

  figure
  hold on
  col = lines(nc);
  legendInfo = cell(1,nc);
  for j = 1:nc
    plot3(Theta{j}(:,col1), Theta{j}(:,col2), Theta{j}(:,col3), '.', 'color', col(j,:));
    xlabel(Labels(col1),'FontSize',16);
    ylabel(Labels(col2),'FontSize',16);
    zlabel(Labels(col3),'FontSize',16);
    grid on;
    title('Evolution of parameters','FontSize',18);
    legendInfo{j} = sprintf('Stage = %d', j-1);
  end
  legend(legendInfo);

else
  error('Not enough input arguments');
end
%% END