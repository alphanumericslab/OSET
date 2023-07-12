function [ARCoefsKF , ARCoefsKS] = time_variant_ar_kalman(data, order, x0, q, R, p0, alpha)
% time_variant_ar_kalman - Time variant auto-regressive (AR) model estimated by Kalman Filter and Kalman Smoother
%
% Syntax: [ARCoefsKF , ARCoefsKS] = time_variant_ar_kalman(data, order, x0, q, R, p0, alpha)
%
% Inputs:
%   data: Template noise used for model training
%   order: AR model order
%   x0: A time-invariant set of AR-coefficients estimated by applying a global AR-model estimation on the entire input signal
%   q: AR coefficients covariance
%   R: Noise variance
%   p0: Covariance of the KF initial state
%   alpha: KF forgetting factor (alpha = 1 for standard KF)
%
% Outputs:
%   ARCoefsKF: AR coefficients estimated by a Kalman Filter
%   ARCoefsKS: AR coefficients estimated by a Kalman Smoother
%
%   Revision History:
%       2006: First release
%       2023: Renamed from deprecated version TimeVariantAR()
% 
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

N = length(data); % Length of the data
A = 1;
Q = q * eye(order);
Wmean = zeros(order, 1);
Vmean = mean(data);
ARCoefsKF = zeros(order, N);
Xminus = x0;
Pminus = p0 * eye(order);

Pbar = zeros(order, order, N);
Phat = zeros(order, order, N);
Xhat = zeros(order, N);
Xbar = zeros(order, N);

% Filtering
for i = 1:N
    Pbar(:, :, i) = Pminus;
    Xbar(:, i) = Xminus;

    if (i < order + 1)
        H = [-data(i - 1:-1:1)' zeros(1, order - i + 1)];
    else
        H = -data(i - 1:-1:i - order)';
    end

    Yminus = H * Xminus + Vmean;
    inov = data(i) - Yminus;

    K = Pminus * H' / (H * Pminus * H' + alpha * R);
    Pplus = ((eye(order) - K * H) * Pminus * (eye(order) - K * H)' + K * R * K') / alpha;
    Xplus = Xminus + K * inov; % A posteriori state estimate

    Xminus = A * Xplus + Wmean; % State update
    Pminus = A * Pplus * A' + Q;

    ARCoefsKF(:, i) = Xplus;
    Phat(:, :, i) = Pplus;
    Xhat(:, i) = Xplus;
end

% Smoothing
PSmoothed = zeros(size(Phat));
X = zeros(size(Xhat));
PSmoothed(:, :, N) = Phat(:, :, N);
X(:, N) = Xhat(:, N);
for k = N - 1:-1:1
    S = Phat(:, :, k) * A' / Pbar(:, :, k + 1);
    X(:, k) = Xhat(:, k) + S * (X(:, k + 1) - Xbar(:, k + 1));
    PSmoothed(:, :, k) = Phat(:, :, k) - S * (Pbar(:, :, k + 1) - PSmoothed(:, :, k + 1)) * S';
end

ARCoefsKS = X;
