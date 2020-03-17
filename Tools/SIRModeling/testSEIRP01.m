% A simulation of the susceptible-exposed-infected-recovered-passed (SEIRP) model for
% epidemic diseases
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
scenario = 'A';

switch scenario
    case 'A'
        T = 4000; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.001*ones(1, K);
        N = 84.0e6;
    case 'B'
        T = 300; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.0*ones(1, K);
        N = 84.0e6;
    case 'C'
        T = 3000; % days
        K = round(T/dt);
        alpha_e = 0.095*linspace(1, 0.1, K);
        alpha_i = 0.005*linspace(1, 0.1, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.001*ones(1, K);
        N = 84.0e6;
    case 'D'
        T = 4000; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.005*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = 0.001*ones(1, K);
        N = 84.0e6;
end

s = zeros(1, K);
e = zeros(1, K);
i = zeros(1, K);
r = zeros(1, K);
d = zeros(1, K);

s(1) = (N - 1)/N;
e(1) = 1/N;
i(1) = 0;
r(1) = 0;
d(1) = 0;

for t = 1 : K - 1
    s(t + 1) = (-alpha_e(t) * s(t) * e(t) - alpha_i(t) * s(t) * i(t) + gamma(t) * r(t)) * dt + s(t);
    e(t + 1) = (alpha_e(t) * s(t) * e(t) + alpha_i(t) * s(t) * i(t) - kappa(t) * e(t) - rho(t) * e(t)) * dt + e(t);
    i(t + 1) = (kappa(t) * e(t) - beta(t) * i(t) - mu(t) * i(t))* dt + i(t);
    r(t + 1) = (beta(t) * i(t) + rho(t) * e(t) - gamma(t) * r(t)) * dt + r(t);
    d(t + 1) = (mu(t) * i(t)) * dt + d(t);
end

t = dt*(0 : K - 1);
figure;
hold on
plot(t, s, 'b', 'linewidth', 3);
plot(t, e, 'c', 'linewidth', 3);
plot(t, i, 'r', 'linewidth', 3);
plot(t, r, 'g', 'linewidth', 3);
plot(t, d, 'k', 'linewidth', 3);
grid
legend('S(t)', 'E(t)', 'I(t)', 'R(t)', 'D(t)');
ylabel('Population Ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');

figure;
hold on
plot(t, i, 'linewidth', 3);
plot(t, e, 'linewidth', 3);
grid
legend('I(t)', 'E(t)');
ylabel('Population ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
% title('The infection rate over time');