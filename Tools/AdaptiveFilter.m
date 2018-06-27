function [e,n] = AdaptiveFilter(x, delay, taps, mu)
% Adaptive Noise-Canceller/Line-Enhancer
%
% Modes of operation:
% Noise Canceller (Signal in periodic noise):
%                                   y = output , n = noise
% Adaptive Line Enhancement Algorithm (Periodic signal in noise):
%                                   y = noise , n = output
% Reza Sameni (C)
%
% Created May 2007
% Modified June 2018

x = x(:);
xref = [zeros(delay+taps-1,1) ; x(1:end-delay)];
w = ones(taps,1);

n = zeros(size(x));
e = zeros(size(x));
for i = 1:length(x)
    if(i > (taps-1))
        xr = xref(i:-1:i-taps+1);
    else
        xr = [xref(i:-1:1) ; zeros(taps-i,1)];
    end
    n(i) = w(:)' * xr;
    e(i) = x(i) - n(i);
    w = w + 2 * mu * e(i) * xr;
end

n = n(:)';
e = e(:)';