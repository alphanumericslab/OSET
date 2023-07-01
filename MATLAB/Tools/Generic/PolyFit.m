function [y, p] = PolyFit(x, fs, N)
% Least-squares polynomial fit over a signal segment
%
% [y, p] = PolyFit(x, fs, N)
%
% Inputs:
%   x: the (noisy) input signal
%   fs: sampling frequency
%   N: polynomial order
% 
% Outputs:
%   y: fitted mode
%   p: polynomial coefficients
%
% The Open-Source Electrophysiological Toolbox (OSET)
% Reza Sameni
% April 2011

x = x(:)';

L = length(x);
n = 0 : L-1;
t = n / fs;

T = zeros(N, L);
TX = zeros(N, 1);
for i = 1 : N
    T(i,:) = t.^(i-1);
    TX(i) = x * T(i,:)';
end
TT = T * T';

p = pinv(TT) * TX;

y = p' * T;