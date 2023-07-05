function y = tanh_saturation(x, ksigma)
% tanh_saturation - Saturates outlier samples using a tanh shape function
%
% Inputs:
%   x: Input data, can be a vector or a matrix (channels x time)
%   ksigma: Scaling factor for the saturation level
%
% Output:
%   y: Saturated data with outliers replaced by the saturation level
%
% Revision History:
%   2020: First release
%   2023: Renamed from deprecated version TanhSaturation()
%
% References:
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

alpha = ksigma * std(x, [], 2); % Compute the scaling factor based on the standard deviation of each channel
alpha = alpha(:); % Convert to a column vector
y = diag(alpha) * tanh(diag(1./alpha) * x); % Scale the input data and apply the tanh function to saturate outliers

% Alternative implementation using broadcasting
% sd = std(x, [], 2);
% alpha = ksigma * sd;
% T = size(x, 2);
% y = alpha(:, ones(1, T)) .* tanh(x ./ (alpha(:, ones(1, T))));

