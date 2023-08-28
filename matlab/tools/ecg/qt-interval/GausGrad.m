function gr = GausGrad(t, p)
% Compute the gradient of a Gaussian function with respect to its parameters.
%
% Inputs:
%   t: Time samples.
%   p: Gaussian parameters [amplitude, width, center].
%
% Output:
%   gr: Gradient of the Gaussian function with respect to its parameters.
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

% Initialize the gradient matrix
gr = zeros(length(p), length(t));
p = p(:); % Ensure p is a column vector

% Iterate through Gaussian parameters
for i = 1 : 3 : length(p)
    % Compute gradient components for amplitude, width, and center
    gr(i, :) = GausVal(t, [1, p(i + 1), p(i + 2)]);
    gr(i + 1, :) = (((t - p(i + 2)).^2) / (p(i + 1)^3)) .* GausVal(t, [p(i), p(i + 1), p(i + 2)]);
    gr(i + 2, :) = ((t - p(i + 2)) / (p(i + 1)^2)) .* GausVal(t, [p(i), p(i + 1), p(i + 2)]);
end
