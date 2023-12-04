function y = bp_filter_dft_windowed(x, fc, bw, wname)
% bp_filter_dft_windowed - Bandpass filter using windowed discrete Fourier transform (DFT) filtering
%
% Usage:
%   y = bp_filter_dft_windowed(x, fc, bw, @wname)
%
% Inputs:
%   x: vector or matrix of input data (channels x samples)
%   fc: normalized center frequency
%   bw: normalized bandwidth
%   @wname: any valid window function name, see help for Matlab's WINDOW function
%
% Output:
%   y: vector or matrix of filtered data (channels x samples)
%
% Note:
%   - fc and bw are the center frequency and bandwidth of the bandpass filter
%     normalized by the sampling frequency
%
%   Revision History:
%       2010: First release
%       2023: Debugged and renamed from deprecated version BPFilter3()
%
%   Reza Sameni, 2010-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

N = size(x, 2); % Number of samples

S = fft(x, N, 2); % Compute the FFT of the input signal

k = floor((fc - bw / 2) * N):ceil((fc + bw / 2) * N); % Indices corresponding to the specified frequency range
k(k < 1) = []; % Remove indices less than 1
k(k > N) = []; % Remove indices greater than N

w = window(wname, length(k)); % Generate the window function
w = w(:)'; % Reshape the window to a row vector
w_matrix = repmat(w, size(S, 1), 1); % Replicate the window row vector for each channel

S(:, k) = w_matrix .* S(:, k); % Apply the window to the specified frequency range
S(:, N - k + 2) = w_matrix .* S(:, N - k + 2); % Apply the window symmetrically for negative frequencies

S1 = zeros(size(S)); % Initialize a matrix for the filtered FFT result
S1(:, k) = S(:, k); % Copy the filtered frequency range
S1(:, N - k + 2) = S(:, N - k + 2); % Copy the filtered negative frequencies symmetrically

y = real(ifft(S1, N, 2)); % Compute the inverse FFT to obtain the filtered signal