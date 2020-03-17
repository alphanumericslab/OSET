% A simulation of the susceptible-infected-recovered (SIR) model for
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

T = 150; % days
dt = 0.1; % simulation time step (in days)
K = round(T/dt);

alpha = 0.5;
beta = 0.05;
gamma = 0.04; % 0.0

N = 84.0e6;
s = zeros(1, K);
i = zeros(1, K);
r = zeros(1, K);

s(1) = (N - 1)/N;
i(1) = 1/N;
r(1) = 0;

for t = 1 : K - 1
   s(t + 1) = (-alpha * s(t) * i(t) + gamma * r(t)) * dt + s(t);
   i(t + 1) = (alpha * s(t) * i(t) - beta * i(t))* dt + i(t);
   r(t + 1) = (beta * i(t) - gamma * r(t)) * dt + r(t);
end

t = dt*(0 : K - 1);
figure;
hold on
plot(t, s, 'b', 'linewidth', 3);
plot(t, i, 'r', 'linewidth', 3);
plot(t, r, 'g', 'linewidth', 3);
grid
legend('S(t)', 'I(t)', 'R(t)');
ylabel('Population Ratio');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');

figure;
hold on
plot(t, N*i, 'linewidth', 3);
grid
ylabel('Population');
xlabel('days');
set(gca, 'fontsize', 16)
set(gca, 'box', 'on');
title('The infection rate over time');
