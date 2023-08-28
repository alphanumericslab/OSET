function x = ecg_gen_from_phase(params, phase)
% ecg_gen_from_phase - Synthetic ECG generator using a Gaussian mixture model.
%
% Syntax:
%   x = ecg_gen_from_phase(params, phase)
%
% Inputs:
%   params: A vector of length 3*L, where L is the number of Gaussian kernels structured as follows:
%               - alpha = params(1:L): Gaussian kernel amplitudes
%               - b = params(L+1:2*L): Gaussian kernel widths (STD)
%               - theta = params(2*L+1:3*L): Gaussian kernel centers
%           Or a structure with fields: params.alpha, params.b and
%               params.theta with similar definitions
%       Params can be generated directly as described above or obtained from
%       functions such as ECGBeatFitter.
%   phase: Cardiac phase signal obtained from phase_calculator or similar functions.
%
% Output:
%   x: Synthetic ECG time-series.
%
% Description:
%   This function generates synthetic ECG signals using a Gaussian mixture model.
%   The signal is generated in a fitted mode, where
%   y(t) = sum_k alpha(k) * exp[-(t - theta_k)^2 / (2 * b_k^2)]
%
% TODO:
%   Return the analytic Jacobians to accelerate the convergence of the nonlinear
%   least squares solver. Added on Feb. 2019
%
% Revision History:
%   2006: First release.
%   2023: Documented and renamed from deprecated version ECGModel
%
% Reza Sameni, 2006-2023
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

% Initialize the synthetic ECG time-series
x = zeros(size(phase));

% Generate the ECG signal using a Gaussian mixture model
for j = 1:length(a)
    % Calculate the phase difference
    dtheta = rem(phase - theta(j) + pi, 2 * pi) - pi;
    
    % Generate Gaussian kernel contribution to the ECG signal
    x = x + alpha(j) .* exp(-dtheta .^ 2 ./ (2 * b(j) .^ 2));
end
