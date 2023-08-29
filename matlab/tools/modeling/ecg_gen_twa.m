function [ecg, phi] = ecg_gen_twa(N, fs, params, phase_bins, varargin)
% ecg_gen_twa - Synthetic T-wave alternans (TWA) ECG generator
%
% Method: Generates synthetic ECG signals by generating Gaussian mixture
%   models with stochastic deviations in the phase domain, followed by
%   scaling the T-wave segment with alternating amplitude fluctuations and
%   transformation from the phase to time domain. The TWA random sequence is
%   controlled by the provided state-transition matrix (STM), which is used
%   to generate an alternating first-order Markov model that switches ECG
%   T-wave amplitudes in the A-B-A-B-... style that is common in TWA,
%   according to the state transition probabilities in the STM. The function
%   generates synthetic ECG signals by first creating a template in the 
%   phase domain using Gaussian mixture modeling with stochastic deviations.
%   It then applies a warping transformation to match the desired ECG length.
%
% Usage:
%   [ecg, phi] = ecg_gen_twa(N, fs, params, phase_bins, warping_order, seed, twa_off)
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
%       params.STM: The state transition matrix, where STM(i,j) represents the
%           probability of transitioning from state i to state j in the next beat.
%           Each row of STM should sum up to 1.
%       params.state0 (optional): initial state of the Markov model for
%           generating TWA. Default is random (1 or 2)
%       params.T_onset: Anticipated T-wave onset in seconds
%       params.T_offset (optional): Anticipated T-wave offset in seconds. Set to 
%           inf if up to the end of the beat is intended. Default: inf
%       params.TWA_shape_factor (optional): T-wave inversion shape factor (determines
%           tanh inversion function slope). Default is 2.0
%   phase_bins: Number of bins in the phase domain for Gaussian template
%   warping_order (optional): Interpolation order for time-warping
%       transformation. Default is 2
%   seed (optional): randomization seed. Default is random (no seed set)
%   twa_off (optional): does not perform TWA. Used for comparison of
%       results with/without TWA. Default is false.
%
% Outputs:
%   ecg: Synthetic ECG signal
%   phi: Corresponding phase vector for the synthetic ECG signal
%
% Revision History:
%   2023: First release.
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if nargin > 4 && ~isempty(varargin{1})
    warping_order = varargin{1};
else
    warping_order = 2;
end

if nargin > 5 && ~isempty(varargin{2})
    seed = varargin{2};
    rng(seed);
end

if nargin > 6 && ~isempty(varargin{3})
    twa_off = varargin{3};
else
    twa_off = false;
end

if ~isfield(params, 'TWA_shape_factor')
    params.TWA_shape_factor = 2.0;
end

if ~isfield(params, 'T_offset')
    params.T_offset = inf;
end

if ~isfield(params, 'state0')
    params.state0 = randi([1, 2]);
end

% Generate phase values
theta = linspace(-pi, pi, phase_bins);
R_peak_phase_index = find(theta >= 0, 1, 'first');

% Initialize variables
signal_len = 0;
n_gmm = length(params.alpha);
ecg = [];
phi = [];
state = params.state0;

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

    % Shape function used to invert the T-wave segment (inverts the polarity of the T-segment)
    p.T_onset = fs * params.T_onset * phase_bins / beat_len;
    p.T_offset = fs * params.T_offset * phase_bins / beat_len;

    if ~twa_off
        shape_window =  tanh_function(phase_bins, ceil(R_peak_phase_index + p.T_onset), params.TWA_shape_factor)...
            - tanh_function(phase_bins, ceil(R_peak_phase_index + p.T_offset), params.TWA_shape_factor);

        state = markov_model_next_state_gen(params.STM, state);
        if state == 1
            twa_high_low_mode = 1;
        else
            twa_high_low_mode = -1;
        end

        shape_window = twa_high_low_mode * shape_window / 2 * params.TWA_fluct_factor + 1;

        % shape_window = ones(1, phase_bins); % For testing purposes only

        ecg_phase_domain_template_shaped = ecg_phase_domain_template .* shape_window;
    else
        ecg_phase_domain_template_shaped = ecg_phase_domain_template;
    end

    % Align the first and last samples
    % find the line between the first and last samples
    baseline = linspace(ecg_phase_domain_template_shaped(1), ecg_phase_domain_template_shaped(end), phase_bins);

    % remove the bias line
    ecg_phase_domain_template_shaped = ecg_phase_domain_template_shaped - baseline;

    % Apply warping transformation to match beat length
    M = warping_transform(phase_bins, beat_len, warping_order);

    % Concatenate the warped template to the ECG and phase matrices
    ecg = cat(1, ecg, M * ecg_phase_domain_template_shaped(:));
    phi = cat(1, phi, M * theta(:));

    % Update the signal length
    signal_len = length(ecg);
end

% Truncate the ECG and phase matrices to the desired length
ecg = ecg(1:N)';
phi = phi(1:N)';
end

function y = tanh_function(len, center, alpha)
t = 1 : len;
y = tanh(alpha * (t - center));
end
