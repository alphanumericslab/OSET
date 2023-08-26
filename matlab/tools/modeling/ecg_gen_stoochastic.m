function [ecg, phi] = ecg_gen_stoochastic(N, fs, varargin)
% ecg_gen_stoochastic - Generate synthetic single-channel ECG with beat-wise stochastic deviations.
%
% Syntax:
%   [ecg, phi] = ecg_gen_stoochastic(N, fs, params)
%   [ecg, phi] = ecg_gen_stoochastic(N, fs, f, f_deviations, alpha, delta_alpha, b, delta_b, theta, delta_theta, theta0)
%
% Inputs:
%   N: Signal length
%   fs: Sampling rate
%   params: A struct containing the following input parameter fields:
%       f: Average heart rate in Hz (BPM/60)
%       f_deviations: Percentage of beat-wise heart rate deviations (Hz)
%       alpha: Amplitudes of Gaussian functions used for ECG modeling
%       delta_alpha: Percentage of amplitude deviations added per beat
%       b: Widths of Gaussian functions used for ECG modeling
%       delta_b: Percentage of Gaussian wave width deviations added per beat
%       theta: Phase of Gaussian functions used for ECG modeling
%       delta_theta: Percentage of Gaussian center deviations added per beat
%       theta0: Initial phase of the synthetic ECG
%   
% Note: An alternative call option is to directly pass the fields of
%       params as inputs (see second form of Syntax above)
% 
% Outputs:
%   ecg: Single-channel synthetic ECG
%   phi: Vector containing the shifted ECG phase
%
% Description:
%   This function generates synthetic single-channel ECG with beat-wise stochastic deviations.
%   The ECG is modeled using Gaussian functions, and the parameters of these functions
%   are modified on a beat-wise basis to introduce variability
%
% References:
%   - Clifford, G. D., Nemati, S., & Sameni, R. (2010). An artificial vector
%     model for generating abnormal electrocardiographic rhythms. Physiological
%     Measurement, 31(5), 595–609.
%   - Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%     Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%     ECG Signals. EURASIP Journal on Advances in Signal Processing, 2007(1).
%   - McSharry, P. E., Clifford, G. D., Tarassenko, L., & Smith, L. A.
%     (2003). A dynamical model for generating synthetic electrocardiogram signals.
%     IEEE Transactions on Biomedical Engineering, 50(3), 289–294.
%
% Revision History:
%   2022: First release.
%   2023: Renamed from the deprecated version SingleChannelECGGeneratorStochastic
%
% Reza Sameni, 2022-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if nargin == 3 && isstruct(varargin{1})
    params = varargin{1};
    f = params.f;
    f_deviations = params.f_deviations;
    alpha = params.alpha;
    delta_alpha = params.delta_alpha;
    b = params.b;
    delta_b = params.delta_b;
    theta = params.theta;
    delta_theta = params.delta_theta;
    theta0 = params.theta0;
elseif nargin == 11
    f = varargin{1};
    f_deviations = varargin{2};
    alpha = varargin{3};
    delta_alpha = varargin{4};
    b = varargin{5};
    delta_b = varargin{6};
    theta = varargin{7};
    delta_theta = varargin{8};
    theta0 = varargin{9};
else
    error('Incorrect number of inputs');
end

w = 2 * pi * f; % Angular frequency
dt = 1 / fs;    % Time step

phi = zeros(1, N); % Initialize phase vector
ecg = zeros(1, N);  % Initialize ECG signal

phi(1) = theta0; % Set initial phase
d_alpha = alpha;     % Initialize amplitude
d_theta = theta;  % Initialize phase
d_b = b;         % Initialize width
n_gmm = length(alpha);
for i = 1 : N - 1
    dtetai = mod(phi(i) - d_theta + pi, 2 * pi) - pi;

    if i == 1
        ecg(i) = sum(d_alpha .* exp(-dtetai .^ 2 ./ (2 * d_b .^ 2)));
    end

    ecg(i + 1) = ecg(i) - dt * sum(w * d_alpha ./ (d_b .^ 2) .* dtetai .* exp(-dtetai .^ 2 ./ (2 * d_b .^ 2)));

    % Next beat
    phi(i + 1) = phi(i) + w * dt;

    if phi(i + 1) > pi % Beat transitions
        phi(i + 1) = phi(i + 1) - 2 * pi;
        d_alpha = alpha .* (1 + (rand(1, n_gmm) - 0.5) * delta_alpha);
        d_theta = theta .* (1 + (rand(1, n_gmm) - 0.5) * delta_theta);
        d_b = b .* max(0, (1 + (rand(1, n_gmm) - 0.5) * delta_b));
        w = 2 * pi * f * max(0, (1 + (rand - 0.5) * f_deviations));
    end
end
