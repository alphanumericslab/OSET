function y = bp_filter_complex_ma_fixed_tap(x, w0, N)
% bp_filter_complex_ma_fixed_tap - Complex frequency domain first-order
%   moving average bandpass filter with a fixed tap
% 
%   y = bp_filter_complex_ma_fixed_tap(x, w0, N)
%   A complex domain moving average bandpass filter.
%
% Inputs:
%   x: input data (matrix or vector)
%   w0: normalized center frequency
%   N: Number of moving average filter taps
%
% Output:
%   y: filtered data
%
% Note:
% - The filter shifts the signal to the left in the frequency domain,
%   applies a moving average filter, and then shifts the signal back to
%   the right in the frequency domain.
%
%   Revision History:
%       2008: First release
%       2023: Renamed from deprecated version BPFilterComplex()
% 
%   Reza Sameni, 2008-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

L1 = size(x, 1);
L2 = size(x, 2);

% Shift the signal to the left in the frequency domain
z = x .* (ones(L1, 1) * exp(-1j * 2 * pi * w0 * (0:L2-1)));

% Filter the signal with an order N moving average (complex-valued)
u = filter(ones(1, N), N, z')';

% Shift the signal back to the right in the frequency domain
y = u .* (ones(L1, 1) * exp(1j * 2 * pi * w0 * (0:L2-1)));
