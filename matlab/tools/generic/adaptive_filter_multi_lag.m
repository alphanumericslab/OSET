function [e, n] = adaptive_filter_multi_lag(x, mindelay, maxdelay, taps, mu)
% adaptive_filter_multi_lag - Adaptive Noise-Canceller/Line-Enhancer with variable delay line length.
%
% Syntax: [e, n] = adaptive_filter_multi_lag(x, mindelay, maxdelay, taps, mu)
%
% Inputs:
%   x: Input signal.
%   mindelay: Minimum delay in samples for the reference signal.
%   maxdelay: Maximum delay in samples for the reference signal.
%   taps: Number of filter taps.
%   mu: Step size or adaptation rate.
%
% Outputs:
%   e: Filter output (enhanced or residual signal).
%   n: Filter output (noise or enhanced signal).
%
% Modes of operation:
%   Noise Canceller (Signal in periodic noise):
%       y = output, n = noise
%   Adaptive Line Enhancement Algorithm (Periodic signal in noise):
%       y = noise, n = output
%
%   Revision History:
%       2007: First release
%       2018: Updated
%       2023: Renamed from deprecated version AdaptiveFilterMultipleDelay()
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET


x = x(:);   % Ensure input is a column vector
delay = mindelay : maxdelay;   % Generate a range of delay values
xref = zeros(length(delay), length(x) + taps - 1);   % Initialize reference signal matrix

% Generate the delayed reference signals for each delay value
for i = 1:length(delay)
    xref(i,:) = [zeros(delay(i) + taps - 1, 1); x(1:end-delay(i))];
end

w = ones(length(delay), taps);   % Initialize filter coefficients matrix

xr = zeros(length(delay), taps);   % Initialize reference signal matrix
n = zeros(length(delay), length(x));   % Initialize noise signal matrix
e = zeros(length(delay), length(x));   % Initialize output signal matrix

% Perform adaptive filtering for each sample of the input signal
for i = 1:length(x)
    for j = 1:length(delay)
        if(i > (taps-1))
            xr(j,:) = xref(j, i:-1:i-taps+1);   % Select the appropriate reference signal
        else
            xr(j,:) = [xref(j, i:-1:1), zeros(1, taps-i)];   % Pad with zeros for initial taps
        end
    end
    n(:, i) = sum(w .* xr, 2);   % Filtered reference signal (noise)
    e(:, i) = x(i) * ones(length(delay), 1) - n(:, i);   % Output signal (residual or enhanced)
    [~, I] = min(abs(e(:, i)));   % Select the index of the minimum absolute value of the output signal
    w = w + 2 * mu * e(I, i) * xr;   % Update filter coefficients using LMS algorithm
end

% Transpose noise signal and output signal to match dimensions
n = n.';
e = e.';

