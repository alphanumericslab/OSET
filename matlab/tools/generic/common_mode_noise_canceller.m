function [x_den, x_den_smoothed, s, alpha, bl] = common_mode_noise_canceller(x, fs, varargin)
%COMMON_MODE_NOISE_CANCELLER Estimates and removes common-mode noise from multichannel data
%
%   [x_den, x_den_smoothed, s, alpha] = common_mode_noise_canceller(x, fs, varargin)
%   estimates and removes common-mode noise from a given multichannel signal
%   matrix `x`, sampled at frequency `fs`. This function implements an iterative
%   algorithm to minimize the squared error between the observed signals and
%   their estimates, produced by common-mode noise removal. The assumption
%   is that a common noise of unknown amplitude has impacted all channels
%   of x.
%
%   Inputs:
%       x - A matrix representing multichannel data, with size N channels x T samples.
%       fs - Sampling frequency of the signal.
%       fc (optional) - Cutoff frequency for low-pass filtering. If not
%              specified, no filtering is applied. If set to 0, the mean
%              of `x` is subtracted per channel.
%       itr (optional) - Number of iterations for the algorithm. Default is 10.
%       lambda (optional) - Regularization parameter for a 2nd order Tikhonov
%              regularization. Default is 100. If set to 0, no
%              regularization is applied.
%       optim_indexes (optional) - Indices of the samples to optimize.
%              Default is all samples.
%       plot_results (optional) - Boolean to enable/disable result plotting.
%              Default is false (no plotting).
%
%   Outputs:
%       x_den - The denoised signal, with common-mode noise removed (N x T).
%       x_den_smoothed - The denoised signal after smoothing with Tikhonov regularization (N x T).
%       s - The estimated common-mode signal (1 x T).
%       alpha - The scaling factor or matrix used in noise cancellation (N x 1).
%       bl - channel-wise average or baseline estimated by high-pass
%           filtering (N x T) or empty, if no `fc` is set 
%
%   Method:
%   This function addresses common-mode noise cancellation in multichannel
%   data by using iterative least squares optimization on a data matrix `X`
%   (N channels by T samples). It aims to minimize the norm of the
%   difference between `X` and its estimate `X_hat` = alpha * s, where
%   `alpha` is an Nx1 vector and `s` is a 1xT vector representing the
%   common-mode signal. The optimization problem is formulated as
%   minimizing the Frobenius norm ||X - alpha * s||^2. The algorithm
%   initializes `s` as the median of `X` across channels and iteratively
%   updates `alpha` and `s` by applying least squares. It allows focusing
%   on specific segments of `X` through `optim_indexes`.
%
%   Example usage:
%       [x_den, x_den_smoothed, s, alpha] = common_mode_noise_canceller(x, fs, 0.5, 15, 200, 10:100, true);
%
%   Revision History:
%       2024: First release
%
%   Reza Sameni, 2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET
%

% Parse optional arguments
if nargin > 2 && ~isempty(varargin{1})
    fc = varargin{1};
else
    fc = [];
end

if nargin > 3 && ~isempty(varargin{2})
    itr = varargin{2};
else
    itr = 10;
end

if nargin > 4 && ~isempty(varargin{3})
    lambda = varargin{3};
else
    lambda = 100;
end

if nargin > 5 && ~isempty(varargin{4})
    optim_indexes = varargin{4};
else
    optim_indexes = 1:size(x, 2);
end

if nargin > 6 && ~isempty(varargin{5})
    plot_results = varargin{5};
else
    plot_results = false;
end

% Baseline subtraction or low-pass filtering, based on `fc`
if ~isempty(fc)
    if fc == 0
        bl = mean(x, 2);
        x = x - bl;
    else
        bl = lp_filter_zero_phase(x, fc/fs);
        x = x - bl;
    end
else
    bl = [];
end

% Initial estimate of common-mode signal `s`
s = median(x, 1);
s = s(:, optim_indexes);

% Iteratively solve for `alpha` and `s`
for k = 1 : itr
    alpha = x(:, optim_indexes) / s;
    s = alpha \ x(:, optim_indexes);
end

% Calculate denoised signal and smoothed version
s = alpha \ x;
x_den = x - alpha * s;

if lambda > 0
    smoothing_order = 2; % 2nd order Tikhonov regularization
    x_den_smoothed = tikhonov_regularization(x_den, smoothing_order, lambda);
else
    x_den_smoothed = x_den;
end

% Optional plotting of results
if plot_results
    figure;
    plot(x');
    hold on;
    plot(s, 'k', 'linewidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Original Signals with Estimated Common-Mode Signal');
    set(gca, 'fontsize', 14)

    figure;
    subplot(211);
    plot(x');
    grid on;
    title('Original Signals');
    xlabel('Time (s)');
    ylabel('Amplitude');
    set(gca, 'fontsize', 14)

    subplot(212);
    plot(x_den');
    grid on;
    title('Denoised Signals');
    xlabel('Time (s)');
    ylabel('Amplitude');
    set(gca, 'fontsize', 14)
end

