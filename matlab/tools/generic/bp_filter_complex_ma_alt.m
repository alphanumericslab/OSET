function y = bp_filter_complex_ma_alt(x,fc,bw,order)
% bp_filter_complex_ma_alt - Bandpass filter using multi-stage moving average filter with
%   frequency shift (an alternative implementation for bp_filter_complex_ma)
%
% Syntax: y = bp_filter_complex_ma_alt(x,fc,bw,order)
%
% Inputs:
%   x: vector or matrix of input data (channels x samples)
%   fc: normalized center frequency (center frequency in Hz / fs)
%   bw: normalized bandwidth (bandwidth in Hz / fs)
%   order: MA filter order
%
% Output:
%   y: vector or matrix of filtered data (channels x samples)
%
% Note:
%   - fc and bw are the center frequency and bandwidth of the bandpass filter
%       normalized by the sampling frequency
%   - The function essentially functions similarly to bp_filter_complex_ma,
%       but has a larger group delay due to using the conv() function in
%       forming the multi-stage MA impulse response 
%
%   Revision History:
%       2010: First release
%       2023: Speedup and renamed from deprecated version BPFilter4()
%
%   Reza Sameni, 2010-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% General parameters
M = size(x,1); % Number of channels
N = size(x,2); % Number of samples
L = round(2/bw); % Length of the moving average filter

% The MA filter impulse response
h0 = ones(1,L)/L;
h = h0;
for i = 1:order-1
    h = conv(h0, h); % Cascade multiple moving average filters
end

% The exponentials used for frequency shifting
n = 0 : N-1; % Sample indices
w = exp(-1j*2*pi*fc*n);
v = exp(1j*2*pi*fc*n);

% Apply frequency shifting and moving average filter to the signal (Real Part)
x1 = x .* w(ones(M,1),:);
y1 = filter(h, 1, x1, [], 2);
z1 = y1 .* v(ones(M,1),:);

% Apply frequency shifting and moving average filter to the signal (Imaginary Part)
x2 = x .* v(ones(M,1),:);
y2 = filter(h, 1, x2, [], 2);
z2 = y2 .* w(ones(M,1),:);

% Combine the real parts of the two filtered signals
y = real(z1 + z2);
