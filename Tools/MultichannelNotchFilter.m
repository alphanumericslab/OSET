function [s, W, A, B, C] = MultichannelNotchFilter(x, fc, Q, fs)
% A multichannel notch filter using frequency domain source separation
% designed for removing powerline and its harmonics from mulctchannel
% recordings
%
% inputs:
% fc: an array of desired notch frequencies
% Q: an array of desired notch Q factors (in forward path)
% fs: the sampling frequency
% 
% Ref: Sameni, Reza, Christian Jutten, and Mohammad B. Shamsollahi. "A 
% deflation procedure for subspace decomposition." IEEE Transactions on 
% Signal Processing 58.4 (2010): 2363-2374.
%
% Reza Sameni
% Feb 2019
%

y = x;
for i = 1:size(x, 1),
    for k = 1:length(fc),
        Wc = fc(k)/(fs/2);
        BW = Wc/Q(k);
        [b_notch_filter, a_notch_filter] = iirnotch(Wc, BW);
        y(i, :) = filtfilt(b_notch_filter, a_notch_filter, y(i, :));
    end
end

B = cov(y');
C = cov(x');

B = (B+B')/2;
C = (C+C')/2;

[V,D] = eig(B, C, 'chol');
[YY, II] = sort(diag(D), 1, 'descend');
W = V(:,II)';
A = pinv(W);
s = W*x;

