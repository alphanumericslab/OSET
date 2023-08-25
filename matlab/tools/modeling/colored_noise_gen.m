function n = colored_noise_gen(sd, Len, fs, beta)
% colored_noise_gen - Generate colored noise by frequency domain filtering of white noise.
%
% Usage:
%   n = colored_noise_gen(sd, Len, fs, beta)
%
% Inputs:
%   sd: Standard deviation of the generated noise.
%   Len: Noise vector length.
%   fs: Sampling rate.
%   beta: Noise color parameter. Assumes the noise spectrum is proportional
%   to 1/f^beta. beta = 0 corresponds to white noise.
%
% Output:
%   n: Colored noise vector
%
% Description:
%   This function generates colored noise by applying frequency domain
%   filtering to white noise. The noise spectrum is modified to follow a
%   1/f^beta shape, creating colored noise.
%
% References:
%   Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%   Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%   ECG Signals. EURASIP Journal on Advances in Signal Processing, 2007(1)
% 
% Revision History:
%   2006: First release
%   2023: Documented and renamed from deprecated version ColoredNoise
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Calculate the length of the Fourier transform
len = Len + ceil(Len / 10);
len = len + (4 - mod(len, 4));

halflen = ceil(len / 2);

% Generate white noise
s = sd * randn(len, 1);

% Compute the power spectrum
S = fft(s, len).^2;

% Frequency vector
f = (0:len-1)' * fs / len;

% Modify the power spectrum based on 1/f^beta relationship
S(2:halflen+1) = S(2:halflen+1) ./ abs(f(2:halflen+1).^beta);
S(halflen+2:len) = S(halflen+2:len) ./ abs(f(halflen:-1:2).^beta);

% Generate the colored noise in the time domain
n = real(ifft(sqrt(S), len, 1));

% Trim the generated noise to the desired length
n = n(halflen - ceil(Len / 2):halflen + floor(Len / 2) - 1);

% Standardize the noise vector
n = sd * (n - mean(n)) / std(n);
n = n(:);
