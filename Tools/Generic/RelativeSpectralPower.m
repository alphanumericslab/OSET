function ratios = RelativeSpectralPower(x, fs, freqs, Q_factor)
% Calculation of relative signal power in narrow frequency bands
% Application: Calculating powerline noise (and harmonics) for total harmonics distortion (THD) analysis 
% Usage:
% ratios = RelativeSpectralPower(x, fs, freqs, Q_factor)
%   Inputs:
%       x: multi-channel input signal (channels x samples)
%       fs: sampling frequency (Hz)
%       freqs: a vector of N frequencies (Hz) (N x 1)
%       Q_factor: the Q-factors of the notch-filters. If Q_factor is a scalar, all frequencies use the same Q_factor 
%   Output:
%       ratios: the relative power of each frequency to the total channel power (channels x N)
%
% Reza Sameni, June 2022
% reza.sameni@gmail.com
% Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET.git

if isscalar(Q_factor)
    Q_factor = repmat(Q_factor, 1, length(freqs));
end

ratios = zeros(size(x, 1), length(freqs));
x_power = sum(x.^2, 2);
for k = 1 : length(freqs)
    W = freqs(k)/(fs/2);
    BW = W/Q_factor(k);
    [b, a] = iirnotch(W, BW);
    x_mains_cancelled = filtfilt(b, a, x')';
    x_mains = x - x_mains_cancelled;
    x_mains_power = sum(x_mains.^2, 2);
    ratios(:, k) = x_mains_power ./ x_power;
end