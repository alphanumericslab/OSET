function y = bp_filter_complex_equiripple(x, fc, bw, factor)
% bp_filter_complex_equiripple - Bandpass filter using an equiripple FIR filter
%   y = bp_filter_complex_equiripple(x, fc, bw, factor)
%   Bandpass filter using an equiripple FIR filter.
%
% Inputs:
%   x: vector or matrix of input data (channels x samples)
%   fc: normalized center frequency
%   bw: normalized bandwidth
%   factor: controls the sharpness of the filter and its convergence (0 < factor < 1)
%
% Output:
%   y: vector or matrix of filtered data (channels x samples)
%
% Note:
% - fc and bw are the center frequency and bandwidth of the bandpass filter
%   normalized by the sampling frequency
%
%
%   Revision History:
%       2014: First release
%       2023: Renamed from deprecated version BPFilter6()
%
%   Reza Sameni, 2014-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% General parameters and dimensions
M = size(x, 1);
N = size(x, 2);
n = 0 : N-1;

% Design the equiripple FIR filter
[nn, fo, mo, w] = firpmord([factor*bw/2 bw/2/factor], [1 0], [0.001 0.1], 1);
h = firpm(nn, fo, mo, w);

% Frequency shifting using complex exponentials
w = exp(-1j*2*pi*fc*n);
v = exp(1j*2*pi*fc*n);

x1 = x .* repmat(w, M, 1); % shift left the signal's center frequency to the origin (baseband)
y1 = filter(h, 1, x1')'; % baseband filter the complex signal
z1 = y1 .* repmat(v, M, 1); % shift the signal to the right back to its original frequency

x2 = x .* repmat(v, M, 1); % shift right the signal's center frequency to the origin (baseband)
y2 = filter(h, 1, x2')'; % baseband filter the complex signal
z2 = y2 .* repmat(w, M, 1); % shift the signal to the left back to its original frequency

% Combine the results and take the real part
y = real(z1 + z2);
