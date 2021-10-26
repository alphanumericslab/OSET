clear
close all
clc

N = 1000000; % noise ensemble length
M = 1000; % number of samples in random sampling batch

sigma = 10.0; % the true noise std
x = sigma * randn(1, N); % the noise ensemble

Itr = 1000; % number of repetitions

n_std_est = zeros(1, Itr);
for k = 1 : Itr
    I = randi(N, [1, M]);
    x_sub = x(I);
    %n_std_est(k) = std(x_sub, 0);
    %     n_std_est(k) = sqrt(var(x_sub));
    mn = mean(x_sub);
    n_std_est(k) = sqrt(sum((x_sub - mn).^2)/(length(x_sub) - 1));
end

figure
hold on
plot(n_std_est)
plot(sigma * ones(1, Itr), 'linewidth', 3);
plot(mean(n_std_est) * ones(1, Itr), 'linewidth', 3);
grid
legend('Random samples', 'True std', 'Mean of random samples');