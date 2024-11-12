function [a_kf_prior, a_kf_post, a_ks, P_kf_prior, P_kf_post, P_ks, kf_gain, eta, a_kf_pred, P_kf_pred, innovations] = tvar_kalman_sys_identify(y, L, varargin)
% tvar_kalman_filter - A Kalamn filter and smoother using time-varying
% autoregressive (TVAR) model.
%
% Syntax: [a_kf_prior, a_kf_post, a_ks, P_kf_prior, P_kf_post, P_ks, kf_gain, eta, a_kf_pred, P_kf_pred, innovations] = tvar_kalman_sys_identify(y, b, a, q, r, gamma, wlen, mode, predict_forward)
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
%   a_kf_post: Vector of denoised signal using the Kalman filter.
%   a_ks: Vector of denoised signal using the Kalman smoother.
%   P_kf_prior: Covariance matrix of the a priori error vector of the Kalman filter.
%   P_kf_post: Covariance matrix of the a posteriori error vector of the Kalman filter.
%   P_ks: Covariance matrix of the a posteriori error vector of the Kalman smoother.
%   kf_gain: Vector of Kalman filter gain.
%   eta: Innovations process running variance normalize by R. Should be arround 1 (used to monitor the KF healthy function).
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

if nargin > 2 && ~isempty(varargin{1})
    Q = varargin{1};
else
    Q = var(y);
end

if nargin > 3 && ~isempty(varargin{2})
    R = varargin{2};
else
    R = var(y);
end

if nargin > 4 && ~isempty(varargin{3})
    gamma = varargin{3};
else
    gamma = 1;
end

if nargin > 5 && ~isempty(varargin{4})
    innov_monitor_win_len = varargin{4};
else
    innov_monitor_win_len = round(0.1 * length(y));
end

% if nargin > 7 && ~isempty(varargin{5})
%     mode = varargin{5};
% else
%     mode = 'obsvr';
% end

if nargin > 6 && ~isempty(varargin{6})
    predict_forward = varargin{6};
else
    predict_forward = 1;
end

% if length(b) >= length(a)
%     error('model is not strictly proper');
% end

%//////////////////////////////////////////////////////////////////////////
% Kalman filter parameters

% Normalize the coefficients
% b = b / a(1);
% a = a / a(1);
% 
% b = b(:);
% b = [zeros(length(a) - 2, 1); b]; % Make the model strictly proper
% a = a(:);
% 
% switch mode
%     case 'controller' % Controller model
%         L = length(a) - 1;
%         A = [-a(2:end)'; eye(L - 1), zeros(L - 1, 1)];
%         B = [1; zeros(L - 1, 1)];
%         H = b';
%     case 'observer' % Observer model
%         L = length(a) - 1;
%         A = [-a(2:end), [eye(L - 1); zeros(1, L - 1)]];
%         B = b;
%         H = [1, zeros(1, L - 1)];
%     otherwise
%         error('Undefined mode');
% end

A = 1.0 * eye(L);
B = 1;

X0 = [-y(1); zeros(L - 1, 1)];
P0 = 1e8 * Q * diag(ones(1, L));%diag([1, .01*ones(1, L - 1)]);

x_pred = X0;
x_prior = X0;
P_prior = P0;
P_pred = P0;
signal_len = length(y);

innovations_running_power = zeros(innov_monitor_win_len, 1) + R;
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

    if n > L
        H = -y(n-1 : -1 : n-L);
    else
        num_zeros_padded = L - n + 1;
        H = [-y(n-1 : -1 : 1), zeros(1, num_zeros_padded)];
    end
    H = H(:)';

    % Store results
    X_prior(:, n) = x_prior;
    P_kf_prior(:, :, n) = P_prior;

    X_pred(:, n) = x_pred;
    P_kf_pred(:, :, n) = P_pred;
    
    % Measurement update (A posteriori updates)
    y_prior = H * x_prior;
    
    K = P_prior * H' / (H * P_prior * H' + R); % Kalman gain
    P_post = (eye(L) - K * H) * P_prior * (eye(L) - K * H)' + K * R * K'; % Stabilized Kalman cov. matrix
    
    innovations(n) = y(n) - y_prior;
    x_post = x_prior + K * innovations(n); % A posteriori state estimate
    
    % Store results
    X_post(:, n) = x_post;
    P_kf_post(:, :, n) = P_post;
    kf_gain(:, n) = K;

    innovations_running_power = [innovations(n) .^ 2; innovations_running_power(1:end - 1)]; % Observation covariance matrix update (for nonstationary signals)
    innovations_var = mean(innovations_running_power);
    eta(n) = innovations_var ./ R;
    R = gamma * R + (1 - gamma) * innovations_var;
    
    % Time update (A priori updates)
    x_prior = A * x_post; % State update
    P_prior = A * P_post * A' + B * Q * B'; % Cov. matrix update

    x_pred = x_prior;
    P_pred = P_prior;
    if predict_forward > 1
        for cntr = 1 : predict_forward - 1
            x_pred = A * x_pred; % State update
            P_pred = A * P_pred * A' + B * Q * B'; % Cov. matrix update
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

a_kf_pred = X_pred(:, :); % Multi-step-ahead prediction results (pre-measurement)
a_kf_prior = X_prior(:, :); % One-step-ahead prediction results (pre-measurement)
a_kf_post = X_post(:, :); % Filtering results (post-measurement)
a_ks = X_smooth(:, :); % Smoothing results (forward-backward)
