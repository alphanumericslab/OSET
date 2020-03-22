% A simulation of multiple scenarios of the susceptible-exposed-infected-recovered-passed (SEIRP) model for
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
scenario = 'B';

switch scenario
    case 'A' % Immunizing disease
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
    case 'B' % Non-immunizing disease
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
    case 'C'
        T = 120; % days
        K = round(T/dt);
        alpha_e = 0.65*linspace(1, 0.01, K);
        alpha_i = 0.005*linspace(1, 0.01, K);
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
    case 'E'
        T = 4000; % days
        K = round(T/dt);
        alpha_e = 0.65*ones(1, K);
        alpha_i = 0.005*ones(1, K);
        kappa = 0.05*ones(1, K);
        rho = 0.08*ones(1, K);
        beta = 0.1*ones(1, K);
        mu = 0.02*ones(1, K);
        gamma = (1/365)*ones(1, K);
        N = 84.0e6;
end

E0 = 1;
[S, E, I, R, P] = SEIRP(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, N, N-E0, E0, 0, 0, 0, T, dt);

% Find eigenvalues of the linearized system (for s(t) close to 1)
A = [alpha_e(1)-kappa(1)-rho(1), alpha_i(1), 0, 0 ; kappa(1), -beta(1)-mu(1), 0, 0 ; rho(1), beta(1), -gamma(1), 0 ; 0, mu(1), 0, 0];
C = [zeros(3, 1) eye(3)]

[V D] = eig(A)
lambda1 = 0;
lambda2 = -gamma(1);
delta = alpha_e(1) - kappa(1) - rho(1);
lambda3 = (delta - beta(1) - mu(1) + sqrt((beta(1) + mu(1) + delta(1))^2 + 4*kappa(1)*alpha_i(1))) / 2;
lambda4 = (delta - beta(1) - mu(1) - sqrt((beta(1) + mu(1) + delta(1))^2 + 4*kappa(1)*alpha_i(1))) / 2;
lambda = [lambda1, lambda2, lambda3, lambda4];
v1 = [0, 0, 0, 1]';
v2 = [0, 0, 1, 0]';
v3 = [1, (lambda3 - delta)/alpha_i(1), (rho(1)*alpha_i(1)+beta(1)*(lambda3-delta))/alpha_i(1)/(lambda3 + gamma(1)), mu(1)*(lambda3 - delta)/lambda3/alpha_i(1)]';
v3 = v3/sqrt(v3'*v3);
v4 = [1, (lambda4 - delta)/alpha_i(1), (rho(1)*alpha_i(1)+beta(1)*(lambda4-delta))/alpha_i(1)/(lambda4 + gamma(1)), mu(1)*(lambda4 - delta)/lambda4/alpha_i(1)]';
v4 = v4/sqrt(v4'*v4);

t = dt*(0 : K - 1);
ii = E0*(lambda3 - delta)*(lambda4 - delta)/(lambda3 - lambda4)*(exp(lambda4*t) - exp(lambda3*t));
ee = E0/(lambda3 - lambda4)*((lambda3 - delta)*exp(lambda4*t) - (lambda4 - delta)*exp(lambda3*t));

figure;
hold on
plot(t, S/N, 'b', 'linewidth', 3);
plot(t, E/N, 'c', 'linewidth', 3);
plot(t, I/N, 'r', 'linewidth', 3);
plot(t, R/N, 'g', 'linewidth', 3);
plot(t, P/N, 'k', 'linewidth', 3);
grid
legend('S(t)', 'E(t)', 'I(t)', 'R(t)', 'P(t)');
ylabel('Population Ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');

figure;
hold on
plot(t, I, 'linewidth', 3);
plot(t, E, 'linewidth', 3);
% plot(t, ii, 'linewidth', 1);
% plot(t, ee, 'linewidth', 1);
grid
legend('I(t)', 'E(t)', 'II(t)', 'EE(t)');
ylabel('Population ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
% title('The infection rate over time');