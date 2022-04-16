function [y, e, P] = RLS(x, d, p, lambda, delta)
%
% Recursive Least Squares Algorithm
% Reference: 
%    Simon S. Haykin, Adaptive Filter Theory, Edition 3, Prentice Hall, 1996
%    (Chapter 13)
%
% Reza Sameni January 2020

x = x(:)';
d = d(:)';
T = length(x);
P = zeros(p, p, T);
K = zeros(p, T);
w = zeros(p, T);
y = zeros(1, T);

P(:, :, 1) = eye(p) / delta;
gamma = 1/lambda;
for n = p : T
    u = x(n : - 1 : n - p + 1)';
    K(:, n) = gamma*P(:, :, n - 1) * u / (1 + gamma*u'*P(:, :, n - 1)*u);
    y(n) = d(n) - w(:, n - 1)' * u;
    w(:, n) = w(:, n - 1) + K(:, n) * y(n);
    P(:, :, n) = gamma * P(:, :, n - 1) - gamma * (K(:, n)*u')*P(:, :, n - 1);
end
e = x - y;
