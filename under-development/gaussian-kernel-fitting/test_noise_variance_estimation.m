clear
close all
clc

N = 10000; % number of noise ensembles
M = 10; % number of samples in each random sampling batch

sigma = 1.0; % the true noise std
x = sigma * randn(1, N); % the noise ensembles with the selected STD

Itr = 100; % number of times to randomly sample from the large random sample set

nvar_est_ML = zeros(1, Itr); % Noise variance ML estimate
nvar_est_UnBiasedML = zeros(1, Itr); % Noise variance ML estimate
nstd_est_BiasCorrected = zeros(1, Itr); % Noise variance estimate with correction term: https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
nvar_est_MAP = zeros(1, Itr); % Noise variance MAP estimate
for k = 1 : Itr
    I = randi(N, [1, M]); % randomly pick M indexes
    x_sub = x(I); % the subset of samples

    % The ML estimate of the subset variance
    mn = mean(x_sub); % the mean of the subset
    nvar_est_UnBiasedML(k) = sum((x_sub - mn).^2)/(M - 1);
    nstd_est_BiasCorrected(k) = sqrt(2 / (M-1)) * gamma(M/2)/gamma((M-1)/2) * sum((x_sub - mn).^2) / M;
    nvar_est_ML(k) = sum((x_sub - mn).^2)/M;

    % The MAP estimate of the subset variance (assuming an exponential distribution of rate lambda = 1 / alpha)
    c = sum(x_sub.^2);
    alpha = 100 * sigma; % alpha should be orders of magnitude larger than the true variance
    nvar_est_MAP(k) = (- M * alpha + sqrt(M^2 * alpha^ 2 + 8 * alpha * c)) / 4.0;
end

nvar_est_UnBiasedML_variance = var(nvar_est_UnBiasedML);
nvar_est_ML_variance = var(nvar_est_ML);
nvar_est_MAP_variance = var(nvar_est_MAP);

nvar_bias_UnBiasedML = abs(sigma^2 - mean(nvar_est_UnBiasedML))/sigma^2;
nvar_bias_ML = abs(sigma^2 - mean(nvar_est_ML))/sigma^2;
nvar_bias_MAP = abs(sigma^2 - mean(nvar_est_MAP))/sigma^2;

disp(['Unbiased ML estimator (bias, variance) = (', num2str(nvar_bias_UnBiasedML), ', ' , num2str(nvar_est_UnBiasedML_variance), ')'])
disp(['Biased ML estimator (bias, variance) = (', num2str(nvar_bias_ML), ', ' , num2str(nvar_est_ML_variance), ')'])
disp(['MAP estimator (bias, variance) = (', num2str(nvar_bias_MAP), ', ' , num2str(nvar_est_MAP_variance), ')'])

figure
hold on
plot(nvar_est_UnBiasedML)
plot(nstd_est_BiasCorrected);
plot(nvar_est_ML)
plot(nvar_est_MAP)
plot(mean(nvar_est_UnBiasedML) * ones(1, Itr), 'linewidth', 3);
plot(mean(nstd_est_BiasCorrected)^2 * ones(1, Itr), 'linewidth', 3);
plot(mean(nvar_est_ML) * ones(1, Itr), 'linewidth', 3);
plot(mean(nvar_est_MAP) * ones(1, Itr), 'linewidth', 3);
plot(sigma^2 * ones(1, Itr), 'linewidth', 3, 'linestyle', '--');
grid
legend('Unbiased ML', 'Bias Corrected', 'ML', 'MAP', 'Unbiased ML Avg', 'Bias Corrected avg', 'ML Avg', 'MAP Avg', 'True VAR');