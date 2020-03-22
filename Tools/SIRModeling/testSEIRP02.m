% A comperative simulation of the susceptible-exposed-infected-recovered-passed (SEIRP) model for
% epidemic diseases under parameter changes
%
% Reza Sameni, March 2020
% reza.sameni@gmail.com
%
% The Open Source Electrophysiological Toolbox, version 3.14, March 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/

clear;
close all;
clc;

dt = 0.1; % simulation time step (in days)
scenario = 'B';

T = 150; % days
K = round(T/dt);
N = 84.0e6;
E0 = 1;

% parameter set 1
alpha_e = 0.6*ones(1, K);
alpha_i = 0.005*ones(1, K);
kappa = 0.05*ones(1, K);
rho = 0.08*ones(1, K);
beta = 0.1*ones(1, K);
mu = 0.02*ones(1, K);
gamma = 0.001*ones(1, K);
[S1, E1, I1, R1, P1] = SEIRP(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, N, N-E0, E0, 0, 0, 0, T, dt);

% parameter set 2
QStartDay = round(30/dt);
QEndDay = round(90/dt);
alpha_e = [0.6*ones(1, QStartDay), 0.1*ones(1, (QEndDay - QStartDay)), 0.4*ones(1, K - QEndDay)];
alpha_i = [0.005*ones(1, QStartDay), 0.001*ones(1, (QEndDay - QStartDay)), 0.001*ones(1, K - QEndDay)];
kappa = 0.05*ones(1, K);
rho = 0.08*ones(1, K);
beta = 0.1*ones(1, K);
mu = 0.02*ones(1, K);
gamma = 0.001*ones(1, K);
[S2, E2, I2, R2, P2] = SEIRP(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, N, N-E0, E0, 0, 0, 0, T, dt);

t = dt*(0 : K - 1);
% figure;
% hold on
% plot(t, S1/N, 'b', 'linewidth', 3);
% plot(t, E1/N, 'c', 'linewidth', 3);
% plot(t, I1/N, 'r', 'linewidth', 3);
% plot(t, R1/N, 'g', 'linewidth', 3);
% plot(t, P1/N, 'k', 'linewidth', 3);
% grid
% legend('S(t)', 'E(t)', 'I(t)', 'R(t)', 'P(t)');
% ylabel('Population Ratio');
% xlabel('days');
% set(gca, 'fontsize', 16)
% set(gca, 'box', 'on');

figure;
subplot(211);
plot(t, I1/N, 'linewidth', 3);
hold on;
plot(t, E1/N, 'linewidth', 3);
plot(t, I2/N, 'linewidth', 3);
plot(t, E2/N, 'linewidth', 3);
grid
legend('I1(t)', 'E1(t)', 'I2(t)', 'E2(t)');
ylabel('Population ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
subplot(212);
plot(t, P1/N, 'linewidth', 3);
hold on;
plot(t, P2/N, 'linewidth', 3);
grid
legend('P1(t)','P2(t)');
ylabel('Population ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
% title('The infection rate over time');