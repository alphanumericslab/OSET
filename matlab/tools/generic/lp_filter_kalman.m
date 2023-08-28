function [y_smoothed, a, y_hat, y_bar] = lp_filter_kalman(y, fc, varargin)
%
% Applies lowpass filtering using a first-order Kalman Smoother.
% 
% Usage:
%   [y_smoothed, a, y_hat, y_bar] = lp_filter_kalman(y, fc, VarWinlen, VarWinlen1, gamma)
%
% Inputs:
%   y: Vector or matrix of input data (channels x samples)
%   fc: -3dB cut-off frequency normalized by the sampling frequency
%   VarWinlen: Length of window used for adapting the variance of the desired signal.
%            Default: 1/10th of the signal length.
%   VarWinlen1: Length of window used for monitoring the whiteness of the innovations signal.
%             Default: 1/10th of the signal length.
%   gamma: Smoothing factor for updating the noise covariance.
%        Default: 1 (no smoothing, identity covariance).
%
% Outputs:
%   y_smoothed: Vector or matrix of filtered data (channels x samples)
%   a: Measure of innovations signal whiteness. 'a' should be close to 1. If
%       'a' oscillates between 0.5 < a < 2, the filter is working well as the
%       innovations signal is relatively white.
%   y_hat: Vector or matrix of Kalman filter posterior esimate of the filtered data (channels x samples)
%   y_bar: Vector or matrix of Kalman filter prior esimate of the filtered data (channels x samples)
% 
% Revision History:
%   2006: First release
%   2023: Renamed from deprecated version KLPFilter
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Parse optional input arguments
if nargin > 2 && ~isempty(varargin{1})
    VarWinlen = varargin{1};
else
    VarWinlen = round(size(y, 2) / 10); % Default: 1/10th of the signal length
end

if nargin > 3 && ~isempty(varargin{2})
    VarWinlen1 = varargin{2};
else
    VarWinlen1 = round(size(y, 2) / 10); % Default: 1/10th of the signal length
end

if nargin > 4 && ~isempty(varargin{3})
    gamma = varargin{3};
else
    gamma = 1; % Default: No smoothing (identity covariance)
end

L = 1;
k = sqrt(2) / 2; % Cut-off value
alpha = (1 - k * cos(2 * pi * fc) - sqrt(2 * k * (1 - cos(2 * pi * fc)) - k^2 * sin(2 * pi * fc)^2)) / (1 - k);
A = alpha;
B = (1 - alpha);
H = 1;
y_smoothed = zeros(size(y));
y_hat = zeros(size(y));
y_bar = zeros(size(y));

Samples = size(y, 2);

a = zeros(size(y));
for i = 1 : size(y, 1)
    ymean = mean(y(i, :));
    y(i, :) = y(i, :) - ymean;

    % Parameters
    bw = filtfilt(1 - alpha, [1, -alpha], y(i, :));
    Q = var(bw);
    R = var(y(i, :) - bw);

    Wmean = 0;
    Vmean = 0;
    X0 = bw(1);
    P0 = 1000 * Q;

    Xminus = X0;
    Pminus = P0;
    Xhat = zeros(1, Samples);
    Phat = zeros(1, Samples);
    Xbar = zeros(1, Samples);
    Pbar = zeros(1, Samples);
    obs = y(i, :);
    mem = zeros(1, VarWinlen) + R;
    mem1 = ones(VarWinlen1, 1);
    %//////////////////////////////////////////////////////////////////////////
    % Filtering
    for k = 1 : Samples
        % Store results
        Xbar(k) = Xminus';
        Pbar(k) = Pminus';

        % Measurement update (A posteriori updates)
        Yminus = H * Xminus + Vmean;

        % Observation available at this sampling time
        K = Pminus * H' / (H * Pminus * H' + R'); % Kalman gain
        Pplus = (eye(L) - K * H) * Pminus * (eye(L) - K * H)' + K * R * K'; % Stabilized Kalman cov. matrix
        innovations = obs(k) - Yminus;
        Xplus = Xminus + K * innovations; % A posteriori state estimate

        Yk = H * Pminus * H' + R;
        mem1 = [innovations.^2 / Yk ; mem1(1:end - 1, :)];
        a(i, k) = mean(mem1);
        mem = [innovations^2 mem(1:end - 1)];
        R = gamma * R + (1 - gamma) * mean(mem);

        % Time update (A priori updates)
        Xminus = A * Xplus + Wmean; % State update
        Pminus = A * Pplus * A' + B * Q * B'; % Cov. matrix update

        % Store results
        Xhat(k) = Xplus';
        Phat(k) = Pplus';
    end

    %//////////////////////////////////////////////////////////////////////////
    % Smoothing
    PSmoothed = zeros(size(Phat));
    X = zeros(size(Xhat));
    PSmoothed(Samples) = Phat(Samples);
    X(Samples) = Xhat(Samples);
    for k = Samples - 1 : -1 : 1
        S = Phat(k) * A' / Pbar(k + 1);
        X(k) = Xhat(k) + S * (X(k + 1) - Xbar(k + 1));
        PSmoothed(k) = Phat(k) - S * (Pbar(k + 1) - PSmoothed(k + 1)) * S';
    end

    y_smoothed(i, :) = X + ymean;
    y_hat(i, :) = Xhat;
    y_bar(i, :) = Xbar;
end

