function [x_kmh] = rk_discrete(diff_eq,x_0,u,h)
%%  [x_kmh] = rk_discrete(diff_eq,x_0,u,h)
%   Continuous to discrete nonlinear system using fourth order Runge-Kutta
%
%   Input data:
%
%   - diff_eq: Handle function with the differential equation
%   - x_0    : State vector at time 't'
%   - u      : Exogenous input at time 't'
%   - h      : Integration step
%
%   Output data:
%
%   - x_kmh  : State vector at time 't+h'
%
%   Source:
%
%   - http://mathworld.wolfram.com/Runge-KuttaMethod.html
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

%% 4th order Runge-Kutta:
k1 = h*diff_eq(x_0,        u);
k2 = h*diff_eq(x_0 + k1/2, u);
k3 = h*diff_eq(x_0 + k2/2, u);
k4 = h*diff_eq(x_0 + k3,   u);

x_kmh = x_0 + (k1 + 2*k2 + 2*k3 + k4)/6;

end
%% END