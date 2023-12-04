function [ecg, phi] = ecg_gen_gmm(phi, theta0, alpha, b, theta)
% ecg_gen_gmm - Generate synthetic ECG using a Gaussian mixture model.
%
% Syntax:
%   [ecg, phi] = ecg_gen_gmm(phi, theta0, alpha, b, theta)
%
% Inputs:
%   phi: ECG phase calculated from real ECG.
%   theta0: Desired phase shift of the synthetic ECG.
%   alpha: Amplitudes of Gaussian functions used for GMM modeling.
%   b: Widths of Gaussian functions used for GMM modeling.
%   theta: Phases of Gaussian functions used for GMM modeling.
%
% Outputs:
%   ecg: Synthetic ECG signal.
%   phi: Vector containing the shifted ECG phase.
%
% Description:
%   This function generates synthetic ECG using a Gaussian mixture model.
%   The ECG is modeled using Gaussian functions with specified amplitudes,
%   widths, and phases.
%
% Revision History:
%   2007: First release.
%   2023: Renamed from the deprecated version SingleChannelECGGenerator.
%
% References:
%   See the references for further details.
%
% Reza Sameni, 2007-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

N = length(phi);

% Shift ECG phase by theta0 and wrap it within [-pi, pi]
phi = mod(phi + theta0 + pi, 2 * pi) - pi;

% Calculate phase differences for Gaussian functions
dtetai = mod(phi(ones(length(theta), 1), :)' - theta(ones(1, N), :) + pi, 2 * pi) - pi;

% Generate synthetic ECG using Gaussian mixture model
ecg = sum(alpha(ones(1, N), :) .* exp(-dtetai .^ 2 ./ (2 * b(ones(1, N), :) .^ 2)), 2)';
