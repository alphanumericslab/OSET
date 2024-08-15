function y = notch_filter_cls(x, ff, gamma)
% notch_filter_cls - A constrained least squares notch filter
%
% Syntax: y = notch_filter_cls(x, ff, gamma)
%
% Inputs:
%   x: input signal (channels x time)
%   ff: mains frequency divided by the sampling frequency
%   gamma: regularization factor
%
% Output:
%   y: the denoised signal (after mains cancellation)
%
% Method:
%   The powerline signal s_k is known to satisfy:
%       I) s_{k+1} + s_{k-1} = 2*cos(2*pi*f_mains/fs) * s_k.
%   Signals corrupted by the powerline noise can be modeled as:
%       II) x_k = s_k + n_k
%   (I) and (II) can be combined to solve the constrained least squares
%   problem per channel: S_opt = argmin(|H*S| + gamma*|X - S|) where H is
%   an oscillator's equation in Toeplitz form (from I), X is the vector form
%   of the samples x_k and S is the vector form of s_k. The norm represents
%   the Frobenius norm. The solution is known to be \hat{S} = (I + gamma*H^T*H)^{-1} * X
%
% NOTE: Due to the large matrix inversion involved, this method is only
%   computationally efficient for relatively short-lengthed signals.
%   Consider using the Kalman notch filter in OSET or a conventional IIR
%   notch filter for long records.
%
% Reference:
%
%   R. Sameni, (2012). A linear Kalman Notch Filter for Power-Line
%       Interference Cancellation. In The 16th CSI International Symposium on
%       Artificial Intelligence and Signal Processing (AISP 2012). 2012 16th
%       CSI International Symposium on Artificial Intelligence and Signal
%       Processing (AISP). IEEE. https://doi.org/10.1109/aisp.2012.6313817
%
%   R. Sameni, Online filtering using piecewise smoothness priors:
%       Application to normal and abnormal electrocardiogram denoising, Signal
%       Processing, Volume 133, 2017, Pages 52-63,
%       https://doi.org/10.1016/j.sigpro.2016.10.019.
%
%   Revision History:
%       2024: First release
%
%   Reza Sameni, 2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

N = size(x, 2); % the signal length in time

H = toeplitz([1 ; zeros(N-3, 1)], [1, -2*cos(2*pi*ff), 1, zeros(1, N-3)]);
% alpha = 1e-3;
% H = toeplitz([1 ; zeros(N-4, 1)], [1, -(alpha + 2*cos(2*pi*ff)), 1 + 2*alpha*cos(2*pi*ff), -alpha, zeros(1, N-4)]);

y = x - ((eye(N) + gamma*(H'*H)) \ x')';

if 1
[U, S, V] = svd(H);
s = diag(S'*S);
end
if 0 % visualize the per-sample time-varient filter response
    A = inv(eye(N) + gamma*(H'*H));
    B = zeros(size(A));
    % A = eye(size(A, 1)) - A;
    for kk = 1 : size(A, 1)
        B(kk, :) = circshift(A(kk, :), size(A, 1)/2 -kk);
    end

    figure
    plot(B');
    grid
end

end
