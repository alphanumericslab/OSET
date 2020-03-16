clear;
close all;
clc;

T = 300; % days
d = 1;%0.00001;

alpha = 0.5/d;
beta = 0.001/d;
gamma = 1e-9/d;

N = 84.0e6;
s = zeros(1, T);
i = zeros(1, T);
r = zeros(1, T);

s(1) = (N - 1)/N;
i(1) = 1/N;
r(1) = 0;

for t = 1 : T - 1
   s(t + 1) = (-alpha * s(t) * i(t) + gamma * r(t)) * d + s(t);
   i(t + 1) = (alpha * s(t) * i(t) - beta * i(t))* d + i(t);
   r(t + 1) = (beta * i(t) - gamma * r(t)) * d + r(t);
end

t = 1 : T;
figure;
hold on
plot(t, N*s, 'b');
plot(t, N*i, 'r');
plot(t, N*r, 'g');
grid
legend('S(t)', 'I(t)', 'R(t)');

figure;
hold on
plot(t, N*i);
grid
legend('I(t)');