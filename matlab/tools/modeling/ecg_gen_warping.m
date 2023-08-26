function [ecg, phi] = ecg_gen_warping(N, fs, params, phase_bins, warping_order)
% ecg_gen_warping - Stochastic synthetic ECG generation using warping transformation from the phase domain.
%
% Method: Generates synthetic ECG signals by generating Gaussian mixture
%   models with stochastic deviations in the phase domain, followed by
%   transformation from the phase to time domain.
% Usage:
%   [ecg, phi] = ecg_gen_warping(N, fs, params, phase_bins, warping_order)
%
% Inputs:
%   N: Desired length of the output ECG signal
%   fs: Sampling rate of the output ECG signal
%   params: A structure containing parameters for ECG generation
%       params.f: Average heart rate in Hz (BPM/60)
%       params.f_deviations: Percentage of beat-wise heart rate deviations (Hz)
%       params.alpha: Amplitudes of Gaussian functions used for ECG modeling
%       params.delta_alpha: Percentage of amplitude deviations added per beat
%       params.b: Widths of Gaussian functions used for ECG modeling
%       params.delta_b: Percentage of Gaussian wave width deviations added per beat
%       params.theta: Phase of Gaussian functions used for ECG modeling
%       params.delta_theta: Percentage of Gaussian center deviations added per beat
%   phase_bins: Number of bins in the phase domain for Gaussian template
%   warping_order: Interpolation order for time-warping transformation
%
% Outputs:
%   ecg: Synthetic ECG signal
%   phi: Corresponding phase vector for the synthetic ECG signal
%
% Description:
%   This function generates synthetic ECG signals by first creating a template
%   in the phase domain using Gaussian mixture modeling with stochastic deviations.
%   It then applies a warping transformation to match the desired ECG length.
%
% Revision History:
%   2023: First release.
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Generate phase values
theta = linspace(-pi, pi, phase_bins);

% Initialize variables
signal_len = 0;
n_gmm = length(params.alpha);
ecg = [];
phi = [];

% Generate the ECG signal iteratively until reaching the desired length
while signal_len < N
    % Determine beat length based on average heart rate and deviations
    beat_len = round(fs / (params.f * max(0, (1 + (rand - 0.5) * params.f_deviations))));

    % Apply stochastic deviations to ECG template parameters
    p.alpha = params.alpha .* (1 + (rand(1, n_gmm) - 0.5) * params.delta_alpha);
    p.theta = params.theta .* (1 + (rand(1, n_gmm) - 0.5) * params.delta_theta);
    p.b = params.b .* max(0, (1 + (rand(1, n_gmm) - 0.5) * params.delta_b));

    % Generate Gaussian template in the phase domain
    ecg_phase_domain_template = gaussian_mixture_model(p, theta);

    % align the first and last samples
    % find the line between x(1) and x(end)
    baseline = linspace(ecg_phase_domain_template(1), ecg_phase_domain_template(end), phase_bins);

    % remove the bias line
    ecg_phase_domain_template = ecg_phase_domain_template - baseline;

    % Apply warping transformation to match beat length
    M = warping_transform(phase_bins, beat_len, warping_order);

    % Concatenate the warped template to the ECG and phase matrices
    ecg = cat(1, ecg, M * ecg_phase_domain_template(:));
    phi = cat(1, phi, M * theta(:));

    % Update the signal length
    signal_len = length(ecg);
end

% Truncate the ECG and phase matrices to the desired length
ecg = ecg(1:N)';
phi = phi(1:N)';
