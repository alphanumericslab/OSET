% Test function for subspace_mapping
%
% Reza Sameni, 2025

close all
clear
clc

% Parameters
ITR = 10;
K = 1000;
dim1 = 32;
dim2 = 3;

% Preallocate
e1 = zeros(1, K);
e2 = zeros(1, K);

tic
% Monte Carlo comparison
for k = 1:K
    A = randn(dim1, dim2);
    B = randn(dim1, dim2);
    [~, ~, e1(k)] = subspace_mapping(A, B, 1, 'rotate');
    [~, ~, e2(k)] = subspace_mapping(A, B, ITR, 'rotate');
end
toc

% Plot results
figure
subplot(2,1,1)
plot(e1, 'bx', 'MarkerSize', 10, 'LineWidth', 1.2)
hold on
plot(e2, 'ro', 'MarkerSize', 10, 'LineWidth', 1.2)
grid on
legend({'Single-shot', sprintf('%d-iteration', ITR)}, 'Location', 'best')
xlabel('Sample')
ylabel('Frobenius Error')
title('Subspace Mapping Convergence')
set(gca, 'FontSize', 13)

subplot(2,1,2)
plot(e1 - e2, 'ko', 'MarkerSize', 8, 'LineWidth', 1.2)
grid on
xlabel('Sample')
ylabel('\Delta Error (Single - Iterative)')
set(gca, 'FontSize', 13)

e1_mean = mean(e1)

e2_mean = mean(e2)
