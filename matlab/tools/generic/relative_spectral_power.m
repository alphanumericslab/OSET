function ratios = relative_spectral_power(x, fs, freqs, Q_factor)
% relative_spectral_power - Calculation of relative signal power in narrow frequency bands
% 
% Syntax: ratios = relative_spectral_power(x, fs, freqs, Q_factor)
%
% Inputs:
%   x: multi-channel input signal (channels x samples)
%   fs: sampling frequency (Hz)
%   freqs: a vector of N frequencies (in Hz) to calculate the relative signal powers (N x 1)
%   Q_factor: the Q-factors of the notch-filters. If Q_factor is a scalar, all frequencies use the same Q_factor 
%
% Output:
%   ratios: the relative power of each frequency to the total channel power (channels x N)
%
% Method: Passes the input signal through notch filters at the given
%   center frequencies and the selected Q-factor. Calculates the ratio
%   of the power of the notch residual over the total signal power
%
% Revision History:
%   2022: First release
%   2023: Renamed from deprecated version RelativeSpectralPower()
%
% Reza Sameni, 2022-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if isscalar(Q_factor)
    Q_factor = repmat(Q_factor, 1, length(freqs));
end

ratios = zeros(size(x, 1), length(freqs));
x_power = sum(x.^2, 2); % Calculate the total power of each channel
for k = 1 : length(freqs) % Loop over the list of frequencies
    W = freqs(k) / (fs / 2); % Normalized frequency
    BW = W / Q_factor(k); % Bandwidth
    [b, a] = iirnotch(W, BW); % Design a notch filter
    x_mains_cancelled = filtfilt(b, a, x')'; % Apply the notch filter to cancel the mains frequency
    x_mains = x - x_mains_cancelled; % Obtain the residual signal after mains cancellation
    x_mains_power = sum(x_mains.^2, 2); % Calculate the power of the residual signal
    ratios(:, k) = x_mains_power ./ x_power; % Calculate the relative power for each channel and frequency
end
