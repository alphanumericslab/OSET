function y = bp_filter_complex_ma(x, fc, bw, order)
% bp_filter_complex_ma applies a zero-phase bandpass filter to the input
%   signal using a multi-stage moving average (MA) filter with complex-valued
%   frequency shifting.
%
% Syntax:
%   y = bp_filter_complex_ma(x, fc, bw, order)
%
% Inputs:
%   x - Input data vector or matrix (channels x samples)
%   fc - Normalized center frequency of the bandpass filter (between 0 and 1)
%   bw - Normalized bandwidth of the bandpass filter (between 0 and 1)
%   order - Order of the MA filter
%
% Output:
%   y - Filtered data vector or matrix (channels x samples)
%
% Note:
% - fc and bw are normalized by the sampling frequency.
% - The filter performs forward-reverse filtering successively. It has
%   zero-phase for even MA filter orders and a phase-lag equal to a single
%   stage MA filter for odd MA filter orders.
% - See alternative implementation: bp_filter_complex_ma_alt
%
%   Revision History:
%       2010: First release
%       2023: Speedup and renamed from deprecated version BPFilter5()
%
%   Reza Sameni, 2010-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% General parameters
M = size(x, 1);     % Number of channels
N = size(x, 2);     % Number of samples

% Simple lowpass moving average filter used as a template
L = round(2 / bw);  % Length of the moving average filter
h = ones(1, L) / L; % Coefficients of the moving average filter

% Exponentials used for frequency shifting
n = 0 : N-1;          % Time index
w = exp(-1j * 2 * pi * fc * n); % Complex exponential for forward filtering
v = exp(1j * 2 * pi * fc * n);  % Complex exponential for reverse filtering

% Real part of the filter
x1 = x .* repmat(w, M, 1); % Apply frequency shifting to the input signal

y1 = x1;
for k = 1 : order
    y1 = filter(h, 1, y1, [], 2); % Apply moving average filter to the signal
    y1 = y1(:, end:-1:1); % Reverse the signal (for reverse filtering)
end
if mod(order, 2) == 1
    y1 = y1(:, end:-1:1); % Reverse the signal again for odd filter order
end
z1 = y1 .* repmat(v, M, 1); % Apply reverse frequency shifting

% Imaginary part of the filter
x2 = x .* repmat(v, M, 1); % Apply reverse frequency shifting to the input signal

y2 = x2;
for k = 1 : order
    y2 = filter(h, 1, y2, [], 2); % Apply moving average filter to the signal
    y2 = y2(:, end:-1:1); % Reverse the signal (for reverse filtering)
end
if mod(order, 2) == 1
    y2 = y2(:, end:-1:1); % Reverse the signal again for odd filter order
end

z2 = y2 .* repmat(w, M, 1); % Apply forward frequency shifting

y = real(z1 + z2); % Combine the real and imaginary parts to obtain the filtered signal
