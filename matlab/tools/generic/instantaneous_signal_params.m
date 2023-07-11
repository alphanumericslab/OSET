function [A, f0, bw, fskew, skew, kurt, f1] = instantaneous_signal_params(signal, fs, wlen, window, th, varargin)
% instantaneous_signal_params - Calculates instantaneous time/frequency domain features of an input signal in causal and non-causal modes
%
% Syntax: [A, f0, bw, fskew, skew, kurt, f1] = instantaneous_signal_params(signal, fs, wlen, window, th, flag)
%
% Inputs:
%   signal: input signal as row vector
%   fs: sampling frequency
%   wlen: sliding window length used for parameter calculation
%   window: a row-wise vector of size wlen used for signal windowing before parameter calculation
%   th: instantaneous amplitude threshold for frequency domain parameter calculation to avoid noisy reports when
%       analytical signal envelope is too weak. Set to 0 to calculate in all amplitudes.
%   flag (optional): Set to 0 for causal and 1 for non-causal calculation. Default = 0 (causal)
%
% Outputs:
%   A: Analytical form signal envelope
%   f0: instantaneous frequency using energy density weighting
%   bw: instantaneous bandwidth using energy density weighting
%   fskew: instantaneous frequency skewness using energy density weighting
%   skew: instantaneous time-domain skewness using energy density weighting
%   kurt: instantaneous time-domain kurtosis using energy density weighting
%   f1: instantaneous frequency using Hilbert transform
%
%   Revision History:
%       2015: First release
%       2020: Added help
%       2023: Renamed from deprecated version InstParams()
% 
%   Reza Sameni, 2015-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if nargin > 5
    flag = varargin{1};
else
    flag = 0; % causal
end

% Make the signals row-wise
signal = signal(:)';
window = window(:)';

A = zeros(size(signal));
f0 = zeros(size(signal));

h = hilbert(signal);
f1 = fs * diff(unwrap(angle(h))) / (2 * pi);
f1 = [f1(1) f1];

bw = zeros(size(signal));
fskew = zeros(size(signal));
skew = zeros(size(signal));
kurt = zeros(size(signal)) + 3; % 3 for Gaussian distribution
N = length(signal); % Signal length

if flag == 0 % Causal
    for i = wlen : N
        x = signal(i - wlen + 1 : i);
        x = window .* x; % Windowing
        NFFT = length(x);
        halfNFFT = round(NFFT / 2);
        X = fft(x, NFFT);
        XX = abs(X(1:halfNFFT)).^2;
        X_energy_density = XX / mean(XX);
        f = (fs / 2) * (0 : halfNFFT - 1) / halfNFFT;
        A(i) = sqrt(mean(abs(x).^2));
        skew(i) = skewness(x);
        kurt(i) = kurtosis(x);
        if A(i) > th
            f0(i) = mean(f .* X_energy_density); % Instantaneous frequency using energy density weighting
            bw(i) = sqrt(mean((f - f0(i)).^2 .* X_energy_density)); % Instantaneous bandwidth using energy density weighting
            fskew(i) = nthroot(mean((f - f0(i)).^3 .* X_energy_density), 3); % Instantaneous frequency skewness using energy density weighting
        end
    end

elseif flag == 1 % Non-causal
    for i = ceil(wlen / 2) : N - floor(wlen / 2)
        x = signal(i - ceil(wlen / 2) + 1 : i + floor(wlen / 2));
        x = window .* x; % Windowing
        NFFT = length(x);
        halfNFFT = round(NFFT / 2);
        X = fft(x, NFFT);
        XX = abs(X(1:halfNFFT)).^2;
        X_energy_density = XX / mean(XX);
        f = (fs / 2) * (0 : halfNFFT - 1) / halfNFFT;
        A(i) = sqrt(mean(abs(x).^2));
        skew(i) = skewness(x);
        kurt(i) = kurtosis(x);
        if A(i) > th
            f0(i) = mean(f .* X_energy_density); % Instantaneous frequency using energy density weighting
            bw(i) = sqrt(mean((f - f0(i)).^2 .* X_energy_density)); % Instantaneous bandwidth using energy density weighting
            fskew(i) = power(mean((f - f0(i)).^3 .* X_energy_density), 1 / 3); % Instantaneous frequency skewness using energy density weighting
        end
    end
end
