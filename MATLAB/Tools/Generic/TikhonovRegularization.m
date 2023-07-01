function y = TikhonovRegularization(x, DiffOrderOrFilterCoefs, lambda)
% Tikhonov Regularization
%
% y = TikhonovRegularization(x, DiffOrderOrFilterCoefs, lambda)
% 
% Solves the constrained least squares problem per channel:
%   y_opt = argmin(|x - y| + lambda*|D*y|)
% Input: x (Channels x Time)
% DiffOrderOrFilterCoefs: Smoothness difference order (if scalar) or
%   smoothing filter coefficients (if vector)
% lambda: regularization (penalty) factor
% 
% Reference: R. Sameni, Online filtering using piecewise smoothness priors:
% Application to normal and abnormal electrocardiogram denoising,
% Signal Processing, Volume 133, 2017, Pages 52-63,
% https://doi.org/10.1016/j.sigpro.2016.10.019.
%
% Reza Sameni, May 2020
% reza.sameni@gmail.com
% Open-Source Electrophysiological Toolbox (OSET)
% https://gitlab.com/rsameni/OSET/


if(length(DiffOrderOrFilterCoefs) == 1)
    h = diff([zeros(1, DiffOrderOrFilterCoefs) 1 zeros(1, DiffOrderOrFilterCoefs)], DiffOrderOrFilterCoefs);
else
    h = DiffOrderOrFilterCoefs;
end

L = length(h);
N = size(x, 2);
D = sparse(toeplitz([h(1) zeros(1, N - L)], [h zeros(1, N - L)]));
F = (lambda*(D'*D) + eye(N));

y = x/F;