function x = gaussian_mixture_model(params, t)
% gaussian_mixture_model - Generate synthetic Gaussian mixture model.
%
% Syntax:
%   x = gaussian_mixture_model(params, t)
%
% Inputs:
%   params: A vector of length 3*L, where L is the number of Gaussian kernels structured as follows:
%               - alpha = params(1:L): Gaussian kernel amplitudes
%               - b = params(L+1:2*L): Gaussian kernel widths (STD)
%               - theta = params(2*L+1:3*L): Gaussian kernel centers
%           Or a structure with fields: params.alpha, params.b and
%               params.theta with similar definitions
% 
%   t: Time/phase vector corresponding to the GMM signal
%
% Outputs:
%   x: Synthetic signal generated using the Gaussian mixture model
%
% Description:
%   This function generates synthetic Gaussian mixture model signals, which
%   can be used for probabolity distribution modeling or ECG modeling with
%   specified amplitudes, widths, and phases.
%
% Revision History:
%   2019: First release
%   2023: Renamed from the deprecated version ECGModelTimeBased
%
% References:
%   See the references for further details.
%
% Reza Sameni, 2019-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if isstruct(params)
    alpha = params.alpha;
    b = params.b;
    theta = params.theta;
else
    % Determine the number of Gaussian kernels
    L = length(params) / 3;

    % Extract parameters
    alpha = params(1:L);
    b = params(L + 1 : 2 * L);
    theta = params(2 * L + 1 : 3 * L);
end

% Initialize the synthetic ECG signal
x = zeros(size(t));

% Generate the synthetic ECG signal using Gaussian kernels
for j = 1 : length(alpha)
    dtheta = t - theta(j);
    x = x + alpha(j) .* exp(-dtheta .^ 2 ./ (2 * b(j) .^ 2));
end
