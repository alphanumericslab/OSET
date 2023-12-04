function y = tikhonov_regularization(x, diff_order_or_filter_coefs, lambda)
% tikhonov_regularization - Tikhonov Regularization.
%
% Syntax: y = tikhonov_regularization(x, diff_order_or_filter_coefs, lambda)
%
% Inputs:
%   x: Input signal (Channels x Time).
%   diff_order_or_filter_coefs: Smoothness difference order (if scalar) or smoothing filter coefficients (if vector).
%   lambda: Regularization (penalty) factor.
%
% Output:
%   y: Regularized signal.
%
% Solves the constrained least squares problem per channel:
%   y_opt = argmin(|x - y| + lambda*|D*y|)
%
% Reference:
%   R. Sameni, Online filtering using piecewise smoothness priors: Application to normal and abnormal electrocardiogram denoising,
%   Signal Processing, Volume 133, 2017, Pages 52-63, https://doi.org/10.1016/j.sigpro.2016.10.019.
%
%   Revision History:
%       2020: First release
%       2023: Renamed from deprecated version TikhonovRegularization
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if (length(diff_order_or_filter_coefs) == 1)
    h = diff([zeros(1, diff_order_or_filter_coefs), 1, zeros(1, diff_order_or_filter_coefs)], diff_order_or_filter_coefs);
else
    h = diff_order_or_filter_coefs;
end

L = length(h);
N = size(x, 2);
D = sparse(toeplitz([h(1), zeros(1, N - L)], [h, zeros(1, N - L)]));   % Create the difference matrix
F = (lambda * (D' * D) + eye(N));   % Form the regularized matrix

y = x / F;   % Solve the regularized least squares problem
