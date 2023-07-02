function y = bp_filter_fft(x, fl, fu)
% bp_filter_fft - Bandpass filter using FFT filtering.
%
% Syntax: y = bp_filter_fft(x, fl, fu)
%
% Inputs:
%   x: Vector or matrix of input data (channels x samples).
%   fl: Normalized lower frequency.
%   fu: Normalized upper frequency.
%
% Output:
%   y: Vector or matrix of filtered data (channels x samples).
%
% Note:
% - fl and fu are the lower and upper frequency ranges of the bandpass filter
%   normalized by the sampling frequency.
% - The filter does not perform any windowing on the data.
%
%   Revision History:
%       2008: First release
%       2023: Renamed from deprecated version BPFilter()
%
%   Reza Sameni, 2008-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

N = size(x, 2);

S = fft(x, N, 2);

k = 1 : ceil(fl * N);
if (~isempty(k))
    S(:, [k, N - k + 2]) = 0;
end

k = floor(fu * N) : ceil(N / 2) + 1;
if (~isempty(k))
    S(:, [k, N - k + 2]) = 0;
end

y = real(ifft(S, N, 2));
