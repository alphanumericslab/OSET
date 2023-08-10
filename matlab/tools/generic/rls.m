function [y, e, P] = rls(x, d, p, lambda, delta)
% rls - Recursive Least Squares Algorithm for adaptive filtering.
%
% Syntax: [y, e, P] = rls(x, d, p, lambda, delta)
%
% Inputs:
%   x: Input signal vector of length T.
%   d: Desired signal vector of length T.
%   p: Order of the filter.
%   lambda: Forgetting factor (0 < lambda <= 1).
%   delta: Regularization factor (delta > 0).
%
% Outputs:
%   y: Output signal vector of length T.
%   e: Error signal vector (d - y) of length T.
%   P: Covariance matrix of filter coefficients, size p x p x T.
%
% Reference: 
%   Simon S. Haykin, Adaptive Filter Theory, Edition 3, Prentice Hall, 1996 (Chapter 13)
%
% Revision History:       
%   2020: First release
%   2023: Renamed from deprecated version RLS()
%
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

x = x(:)';      % Convert x to a row vector
d = d(:)';      % Convert d to a row vector
T = length(x);  % Length of input signal
P = zeros(p, p, T); % Covariance matrix of filter coefficients
K = zeros(p, T);    % Kalman gain matrix
w = zeros(p, T);    % Filter coefficients matrix
y = zeros(1, T);    % Output signal vector

P(:, :, 1) = eye(p) / delta;    % Initialize the covariance matrix at time index 1
gamma = 1 / lambda;             % Calculate the gamma parameter based on the forgetting factor
for n = p : T
    u = x(n : - 1 : n - p + 1)';    % Generate the input vector
    K(:, n) = gamma * P(:, :, n - 1) * u / (1 + gamma * u' * P(:, :, n - 1) * u);    % Update the Kalman gain
    y(n) = d(n) - w(:, n - 1)' * u;    % Calculate the output signal
    w(:, n) = w(:, n - 1) + K(:, n) * y(n);    % Update the filter coefficients
    P(:, :, n) = gamma * P(:, :, n - 1) - gamma * (K(:, n) * u') * P(:, :, n - 1);    % Update the covariance matrix
end

e = d - y;  % Calculate the error signal


