% This script tests the goodness of fit of gaussian functions on an ECG

clear;
close all;
clc;


% load('SampleECG1kHz1.mat');
load('SampleECG1kHz2.mat');
data = data';
fs = 1000;
ch = 2;
ksigma = 4;

% x = LPFilter(data - LPFilter(data, 3/fs), 250/fs);
x = data(ch, :) - LPFilter(data(ch, :), 2/fs);
y = -[0, diff(x)]./x;

z_poly = zeros(size(y));
err_poly = zeros(size(y));

z_tikhonov = zeros(size(y));
err_tikhonov = zeros(size(y));

poly_order = 1;
wlen = 5;
if(mod(wlen, 2) == 0)
   error('wlen must be odd'); 
end
wlen_half = floor(wlen);
%     err(k) = norm(yy - yy_poly_fit);
for k = wlen_half + 1 : length(y)-wlen_half
    yy = y(k - wlen_half : k + wlen_half);
    [yy_poly_fit, p] = PolyFit(yy, fs, poly_order);
    err_poly(k) = median(abs((yy - yy_poly_fit)));
    z_poly(k) = yy_poly_fit(wlen_half + 1);

    yy_tikhonov = TikhonovRegularization(yy, [1 -1], 0.1);
    err_tikhonov(k) = median(abs((yy - yy_tikhonov)));
    z_tikhonov(k) = yy_tikhonov(wlen_half + 1);
end

figure
plot(x);
grid

figure
plot(x(1:end-1), -diff(x), 'b.');
grid
xlabel('x');
ylabel('dx/dt');

lgnd = {};
figure
hold on
plot(x); lgnd = cat(1, lgnd, 'x');
plot(y); lgnd = cat(1, lgnd, 'gaussian slope');
plot(z_poly); lgnd = cat(1, lgnd, 'poly fit on gaussian slope');
plot(z_tikhonov); lgnd = cat(1, lgnd, 'tikhonov reg. fit on gaussian slope');
% plot(err, 'm'); lgnd = cat(1, lgnd, 'gaussian slope poly fit error');
grid
xlabel('t');
legend(lgnd);
set(gca, 'fontsize', 16);

figure
hold on
plot(TanhSaturation(y, ksigma), 'k');
plot(TanhSaturation(y, ksigma), 'r.');
plot(x, 'b');
grid
xlabel('t');