function [y_kf_prior, y_kf_post, y_ks, P_kf_prior, P_kf_post, P_ks, kf_gain, eta, y_kf_pred, P_kf_pred, r] = arma_kalman_filter(y, b, a, varargin)
% arma_kalman_filter - Removing lowpass noise using a first order linear Kalman filter and smoother.
%
% Syntax: [y_kf_prior, y_kf_post, y_ks, P_kf_prior, P_kf_post, P_ks, kf_gain, eta, y_kf_pred, P_kf_pred] = arma_kalman_filter(y, b, a, q, r, gamma, wlen, mode, predict_forward)
%
% Inputs:
%   y: Vector of noisy signals contaminated with lowpass noise.
%   b: Numerator of the model transfer function.
%   a: Denominator of the model transfer function.
%   q: Covariance of the model errors. By default: var(y).
%   r: Covariance of the power-line noise. By default: var(y).
%   gamma: Nonstationarity adaptation parameter. 0 <= gamma <= 1 and by default gamma = 1.
%   wlen: Nonstationarity adaptation window length. By default wlen = length(y)/10.
%   mode: Filtering mode ('obsvr' for observer model, 'contr' for controller model). By default 'obsvr'.
%
% Outputs:
%   y_kf_post: Vector of denoised signal using the Kalman filter.
%   y_ks: Vector of denoised signal using the Kalman smoother.
%   P_kf_prior: Covariance matrix of the a priori error vector of the Kalman filter.
%   P_kf_post: Covariance matrix of the a posteriori error vector of the Kalman filter.
%   P_ks: Covariance matrix of the a posteriori error vector of the Kalman smoother.
%   kf_gain: Vector of Kalman filter gain.
%   eta: Innovations process running variance normalize by r. Should be arround 1 (used to monitor the KF healthy function).
%
% Note: The ARMA model is assumed to be in the following form: y(z)/w(z) =
% B(z)/A(Z), where: B(z) = b1 + b2*z^{-1} + ... + bm*z^{-(m-1)} and A(z) = a1 + a2*z^{-1} + ... + an*z^{-(n-1)}
%
%   Revision History:
%       2008: First release
%       2023: Renamed from deprecated version KalmanARFilter()
%       2024: Updated help
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

%//////////////////////////////////////////////////////////////////////////

signal_len = length(y); % input signal length

% Input arguments
if nargin > 3 && ~isempty(varargin{1})
    q = varargin{1};
else
    q = var(y);
end

if nargin > 4 && ~isempty(varargin{2})
    r = varargin{2};
else
    r = var(y);
end

if nargin > 5 && ~isempty(varargin{3})
    gamma = varargin{3};
else
    gamma = 1;
end

if nargin > 6 && ~isempty(varargin{4})
    innov_monitor_win_len = varargin{4};
else
    innov_monitor_win_len = ceil(0.01 * signal_len);
end

if nargin > 7 && ~isempty(varargin{5})
    mode = varargin{5};
else
    mode = 'observer';
end

if nargin > 8 && ~isempty(varargin{6})
    predict_forward = varargin{6};
else
    predict_forward = 1;
end

if length(b) >= length(a)
    error('model is not strictly proper');
end

if isscalar(r)
    r = repmat(r, 1, signal_len);
else
    if length(r) ~= signal_len
        error('r should either be a scalar or a vector with the same length as the input signal');
    end

    if gamma ~= 1
        % gamma = 1;
        % warming('When r is a vector, no noise variance adaptation is performed. Resetting gamma to 1.');
        warning('r is a vector and gamma ~=1; r will be adapted over time.');
    end
end
r = r(:)';

%//////////////////////////////////////////////////////////////////////////
% Kalman filter parameters

% Normalize the coefficients
b = b / a(1);
a = a / a(1);

b = b(:);
b = [zeros(length(a)-2, 1); b]; % Make the model strictly proper
a = a(:);

L = length(a) - 1;
switch mode % See Page 50: T. Kailath, Linear Systems, 1980
    case 'corrected' % Corrected
        A = [-a(2:end)';...
            eye(L - 1), zeros(L - 1, 1)];
        B = [1; zeros(L - 1, 1)];
        H = -a(2:end)';
    case 'controller' % Controller model
        A = [-a(2:end)';...
            eye(L - 1), zeros(L - 1, 1)];
        B = [1; zeros(L - 1, 1)];
        H = b';
    case 'observer' % Observer model
        A = [-a(2:end), [eye(L - 1); zeros(1, L - 1)]];
        B = b;
        H = [1, zeros(1, L - 1)];
    otherwise
        error('Undefined mode');
end

X0 = [median(y, "omitnan") ; zeros(L-1, 1)];
P0 = diag([1e2*q, q*ones(1, L - 1)]);

x_pred = X0;
x_prior = X0;
P_prior = P0;
P_pred = P0;

innovations_running_power = r(1 : innov_monitor_win_len)';
innovations = zeros(1, signal_len);

X_pred = zeros(L, signal_len);
P_kf_pred = zeros(L, L, signal_len);

X_post = zeros(L, signal_len);
P_kf_post = zeros(L, L, signal_len);

X_prior = zeros(L, signal_len);
P_kf_prior = zeros(L, L, signal_len);

kf_gain = zeros(L, signal_len);
eta = zeros(1, signal_len);

%//////////////////////////////////////////////////////////////////////////
% Forward Filtering Stage
for n = 1 : signal_len
    % Store results
    X_prior(:, n) = x_prior;
    P_kf_prior(:, :, n) = P_prior;

    X_pred(:, n) = x_pred;
    P_kf_pred(:, :, n) = P_pred;

    % Measurement update (A posteriori updates)
    y_prior = H * x_prior;

    if ~isnan(y(n)) && ~isinf(r(n))% normal update
        innovations(n) = y(n) - y_prior;
        innovations_var_nominal = H * P_prior * H' + r(n);

        K = P_prior * H' / innovations_var_nominal; % Kalman gain

        x_post = x_prior + K * innovations(n); % A posteriori state estimate
        P_post = (eye(L) - K * H) * P_prior * (eye(L) - K * H)' + K * r(n) * K'; % Stabilized Kalman cov. matrix
    else
        innovations(n) = nan;
        innovations_var_nominal = inf;
        K = zeros(L, 1);
        x_post = x_prior;
        P_post = P_prior;
    end

    % Store results
    X_post(:, n) = x_post;
    P_kf_post(:, :, n) = P_post;
    kf_gain(:, n) = K;

    innovations_running_power = [innovations(n) .^ 2; innovations_running_power(1:end - 1)]; % Observation covariance matrix update (for nonstationary signals)
    innovations_var_practical = mean(innovations_running_power, "omitnan");
    eta(n) = innovations_var_practical / innovations_var_nominal;
    if isfinite(innovations_var_practical) && n < signal_len && gamma ~= 1
        r(n+1) = gamma * r(n+1) + (1 - gamma) * max(innovations_var_practical - H * P_prior * H', 0);
    end

    % Time update (A priori updates)
    x_prior = A * x_post; % State update
    P_prior = A * P_post * A' + B * q * B'; % Cov. matrix update

    x_pred = x_prior;
    P_pred = P_prior;
    if predict_forward > 1
        % disp('Predicting:');
        for cntr = 1 : predict_forward - 1
            % x_pred';
            x_pred = A * x_pred; % State update
            P_pred = A * P_pred * A' + B * q * B'; % Cov. matrix update
        end
    end

end

%//////////////////////////////////////////////////////////////////////////
% Backward Smoothing Stage
P_ks = zeros(size(P_kf_post));
P_ks(:, :, signal_len) = P_kf_post(:, :, signal_len);
X_smooth = zeros(size(X_post));
X_smooth(:, signal_len) = X_post(:, signal_len);
for n = signal_len - 1 : -1 : 1
    S = P_kf_post(:, :, n) * A' / P_kf_prior(:, :, n + 1);
    X_smooth(:, n) = X_post(:, n) + S * (X_smooth(:, n + 1) - X_prior(:, n + 1));
    P_ks(:, :, n) = P_kf_post(:, :, n) - S * (P_kf_prior(:, :, n + 1) - P_ks(:, :, n + 1)) * S';
end

y_kf_pred = X_pred(1, :); % Multi-step-ahead prediction results (pre-measurement)
y_kf_prior = X_prior(1, :); % One-step-ahead prediction results (pre-measurement)
y_kf_post = X_post(1, :); % Filtering results (post-measurement)
y_ks = X_smooth(1, :); % Smoothing results (forward-backward)
