function [e,n] = AdaptiveFilterMultipleDelay(x, mindelay, maxdelay, taps, mu)
% Adaptive Noise-Canceller/Line-Enhancer with variable delay line length
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
delay = mindelay : maxdelay;
xref = zeros(length(delay), length(x) + taps - 1);
for i = 1:length(delay)
    xref(i,:) = [zeros(delay(i)+taps-1,1) ; x(1:end-delay(i))];
end
w = ones(length(delay),taps);

xr = zeros(length(delay), taps);
n = zeros(length(delay), length(x));
e = zeros(length(delay), length(x));
for i = 1 : length(x)
    for j = 1:length(delay)
        if(i > (taps-1))
            xr(j,:) = xref(j, i:-1:i-taps+1);
        else
            xr(j,:) = [xref(j, i:-1:1) , zeros(1, taps-i)];
        end
    end
    n(:, i) = sum(w .* xr, 2);
    e(:, i) = x(i)*ones(length(delay), 1) - n(:, i);
    I = find(abs(e(:, i)) == min(abs(e(:, i))), 1);
    w = w + 2 * mu * e(I, i) * xr;
end
