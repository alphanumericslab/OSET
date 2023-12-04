clear
close all
clc


DiffOrder = 2;
wlen = 100;% in samples
lambda = 200;

h = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder);
L = length(h);
D = toeplitz([h(1) zeros(1, wlen - L)], [h zeros(1, wlen - L)]);

figure
mesh(pinv(eye(size(D,2)) + (lambda)^2*(D'*D)));
% mesh(pinv(eye(size(Dd,2)) + (lambda)^2*(Dd'*Dd)));
set(gca, 'fontsize', 16);
axis tight

