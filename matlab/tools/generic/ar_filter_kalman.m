function [y1, y2, Pbar, Phat, PSmoothed, Kgain] = ar_filter_kalman(x, b, a, varargin)
% ar_filter_kalman - Removing lowpass noise using a first order linear Kalman filter and smoother.
%
% Syntax: [y1, y2, Pbar, Phat, PSmoothed, Kgain] = ar_filter_kalman(x, b, a, q, r, gamma, wlen, mode)
%
% Inputs:
%   x: Vector of noisy signals contaminated with lowpass noise.
%   b: Numerator of the model transfer function.
%   a: Denominator of the model transfer function.
%   q: Covariance of the model errors. By default: var(x).
%   r: Covariance of the power-line noise. By default: var(x).
%   gamma: Nonstationarity adaptation parameter. 0 <= gamma <= 1 and by default gamma = 1.
%   wlen: Nonstationarity adaptation window length. By default wlen = length(x)/10.
%   mode: Filtering mode ('obsvr' for observer model, 'contr' for controller model). By default 'obsvr'.
%
% Outputs:
%   y1: Vector of denoised signal using the Kalman filter.
%   y2: Vector of denoised signal using the Kalman smoother.
%   Pbar: Covariance matrix of the a priori error vector of the Kalman filter.
%   Phat: Covariance matrix of the a posteriori error vector of the Kalman filter.
%   PSmoothed: Covariance matrix of the a posteriori error vector of the Kalman smoother.
%   Kgain: Vector of Kalman filter gain.
%
%   Revision History:
%       2008: First release
%       2023: Renamed from deprecated version KalmanARFilter()
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

%//////////////////////////////////////////////////////////////////////////
% Input arguments

if (nargin > 3 && ~isempty(varargin{1}))
    Q = varargin{1};
else
    Q = var(x);
end

if (nargin > 4 && ~isempty(varargin{2}))
    R = varargin{2};
else
    R = var(x);
end

if (nargin > 5 && ~isempty(varargin{3}))
    gamma = varargin{3};
else
    gamma = 1;
end

if (nargin > 6 && ~isempty(varargin{4}))
    VarWinlen = varargin{4};
else
    VarWinlen = round(length(x) / 10);
end

if (nargin > 7 && ~isempty(varargin{5}))
    mode = varargin{5};
else
    mode = 'obsvr';
end

if (length(b) >= length(a))
    error('model is not strictly proper');
end

%//////////////////////////////////////////////////////////////////////////
% Kalman filter parameters

% Normalize the coefficients
b = b / a(1);
a = a / a(1);

b = b(:);
b = [zeros(length(a) - 2, 1); b]; % Make the model strictly proper
a = a(:);

if (strcmp(mode, 'contr')) % Controller model
    L = length(a) - 1;
    A = [-a(2:end)'; eye(L - 1) zeros(L - 1, 1)];
    B = [1; zeros(L - 1, 1)];
    H = b';
else % Observer model
    L = length(a) - 1;
    A = [-a(2:end) [eye(L - 1); zeros(1, L - 1)]];
    B = b;
    H = [1 zeros(1, L - 1)];
end

X0 = [x(1); zeros(L - 1, 1)];
P0 = 1000 * Q * diag([1 .01 * ones(1, L - 1)]);

Xminus = X0;
Pminus = P0;
Samples = length(x);

mem2 = zeros(VarWinlen, 1) + R;
Xhat = zeros(L, Samples);
innovations = zeros(1, Samples);
Phat = zeros(L, L, Samples);
Xbar = zeros(L, Samples);
Pbar = zeros(L, L, Samples);
Kgain = zeros(L, Samples);

%//////////////////////////////////////////////////////////////////////////
% Forward Filtering Stage
for k = 1 : Samples
    % Store results
    Xbar(:, k) = Xminus;
    Pbar(:, :, k) = Pminus;
    
    % Measurement update (A posteriori updates)
    Yminus = H * Xminus;
    
    K = Pminus * H' / (H * Pminus * H' + R); % Kalman gain
    Pplus = (eye(L) - K * H) * Pminus * (eye(L) - K * H)' + K * R * K'; % Stabilized Kalman cov. matrix
    
    innovations(k) = x(k) - Yminus;
    Xplus = Xminus + K * innovations(k); % A posteriori state estimate
    
    mem2 = [innovations(k) .^ 2; mem2(1:end - 1)]; % Observation covariance matrix update (for nonstationary signals)
    R = gamma * R + (1 - gamma) * mean(mem2);
    
    % Time update (A priori updates)
    Xminus = A * Xplus; % State update
    Pminus = A * Pplus * A' + B * Q * B'; % Cov. matrix update
    
    % Store results
    Xhat(:, k) = Xplus;
    Phat(:, :, k) = Pplus;
    Kgain(:, k) = K;
end

%//////////////////////////////////////////////////////////////////////////
% Backward Smoothing Stage
PSmoothed = zeros(size(Phat));
X = zeros(size(Xhat));
PSmoothed(:, :, Samples) = Phat(:, :, Samples);
X(:, Samples) = Xhat(:, Samples);
for k = Samples - 1 : -1 : 1
    S = Phat(:, :, k) * A' / Pbar(:, :, k + 1);
    X(:, k) = Xhat(:, k) + S * (X(:, k + 1) - Xbar(:, k + 1));
    PSmoothed(:, :, k) = Phat(:, :, k) - S * (Pbar(:, :, k + 1) - PSmoothed(:, :, k + 1)) * S';
end

y1 = Xhat(1, :); % Filtering results
y2 = X(1, :); % Smoothing results
