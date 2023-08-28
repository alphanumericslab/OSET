function y = bp_filter_dft(x, fl, fu)
% bp_filter_dft - Bandpass filter using discrete Fourier transform (DFT) filtering.
%
% Syntax: y = bp_filter_dft(x, fl, fu)
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
%       2023: Debugged and renamed from deprecated version BPFilter
%
%   Reza Sameni, 2008-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if fl > 0.5
    error('fl may not exceed 0.5');
end
if fu > 0.5
    error('fu may not exceed 0.5');
end
if fl > fu
    error('fl may not exceed fu');
end


N = size(x, 2);
S = fft(x, N, 2);

k_cut_low = floor(fl * N);
S(:, [1 : k_cut_low, (N - k_cut_low + 2) : N]) = 0;

k_cut_high = ceil(fu * N);
if k_cut_high > 0
    S(:, k_cut_high + 1 : N - k_cut_high + 1) = 0;
else
    S = 0;
end

y = real(ifft(S, N, 2));
