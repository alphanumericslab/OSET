function [e, n] = adaptive_filter(x, delay, taps, mu)
% adaptive_filter - Adaptive Noise-Canceller/Line-Enhancer.
%
% Syntax: [e, n] = adaptive_filter(x, delay, taps, mu)
%
% Inputs:
%   x: Input signal.
%   delay: Delay in samples for the reference signal.
%   taps: Number of filter taps.
%   mu: adaptation rate.
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
%       2023: Renamed from deprecated version AdaptiveFilter()
%
%   Reza Sameni, 2007-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET


x = x(:);   % Ensure input is a column vector
xref = [zeros(delay+taps-1, 1); x(1:end-delay)];   % Delayed reference signal
w = ones(taps, 1);   % Initialize filter coefficients

n = zeros(size(x));   % Initialize noise signal
e = zeros(size(x));   % Initialize output signal

for i = 1:length(x)
    if(i > (taps-1))
        xr = xref(i:-1:i-taps+1);   % Select the appropriate reference signal
    else
        xr = [xref(i:-1:1); zeros(taps-i, 1)];   % Pad with zeros for initial taps
    end
    n(i) = w(:)' * xr;   % Filtered reference signal (noise)
    e(i) = x(i) - n(i);   % Output signal (residual or enhanced)
    w = w + 2 * mu * e(i) * xr;   % Update filter coefficients using LMS algorithm
end

n = n(:)';   % Transpose noise signal to match dimensions
e = e(:)';   % Transpose output signal to match dimensions
