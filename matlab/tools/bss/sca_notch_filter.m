function [s, W, A, B, C] = sca_notch_filter(x, fc, Q, fs)
%
% Multichannel notch filter using Spectral Component Analysis (SCA). Designed
% for removing powerline and its harmonics from multichannel recordings.
%
% Usage:
%   [s, W, A, B, C] = sca_notch_filter(x, fc, Q, fs);
% 
% Inputs:
%   x: Multichannel input signal (each row represents a channel)
%   fc: Array of desired notch frequencies
%   Q: Array of desired notch Q factors (in forward path)
%   fs: Sampling frequency
%
% Outputs:
%   s: Separated signal after notch filtering
%   W: Mixing matrix learned by SCA
%   A: Unmixing matrix (inverse of W)
%   B: Covariance matrix of the filtered signal
%   C: Covariance matrix of the input signal
%
% Reference:
%   Sameni, Reza, Christian Jutten, and Mohammad B. Shamsollahi.
%   "A deflation procedure for subspace decomposition." IEEE Transactions on
%   Signal Processing 58.4 (2010): 2363-2374.
%
% Revision History:
%   2019: First release
%   2023: Removed channel loop and used matrix form; renamed from deprecated version MultichannelNotchFilter
%
% Reza Sameni, 2019-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Initialize the output variable
y = x;

% Loop through each notch frequency
for k = 1 : length(fc)
    % Normalized cutoff frequency and bandwidth
    Wc = fc(k)/(fs/2);
    BW = Wc/Q(k);

    % Design the notch filter
    [b_notch_filter, a_notch_filter] = iirnotch(Wc, BW);

    % Apply the notch filter to all channels
    y = filtfilt(b_notch_filter, a_notch_filter, y')';
end

% Calculate the covariance matrices of filtered and input signals
B = cov(y');
C = cov(x');

% Ensure the covariance matrices are symmetric
B = (B + B') / 2;
C = (C + C') / 2;

% Perform subspace decomposition using eigenvalue decomposition
[V, D] = eig(B, C, 'chol');

% Sort eigenvalues in descending order and rearrange eigenvectors accordingly
[~, II] = sort(diag(D), 1, 'descend');
W = V(:, II)';
A = pinv(W);

% Separate the sources using unmixing matrix
s = W * x;
