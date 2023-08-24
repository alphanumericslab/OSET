function g = GausVal(t, p)
% Compute the Gaussian function values.
%
% Inputs:
% t: Time samples.
% p: Gaussian parameters [amplitude, width, center]. Can be a matrix where each
%    column corresponds to a set of parameters.
%
% Output:
% g: Values of the Gaussian function evaluated at time samples 't'.
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%   Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions
%   on Signal Processing 70 (2022): 3181-3192.
%
% Revision History:
%   2021: First release
%
% Davood Fattahi (fattahi.d@gmail.com), 2021
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

t = t(:); % Ensure 't' is a column vector

if isvector(p)
    p = reshape(p, 3, []); % Reshape 'p' into a 3xN matrix if it's a vector
end

% Compute the Gaussian function values
g = sum(p(1, :) .* exp(-((repmat(t, 1, size(p, 2)) - p(3, :)).^2) ./ (2 * p(2, :).^2)), 2);
