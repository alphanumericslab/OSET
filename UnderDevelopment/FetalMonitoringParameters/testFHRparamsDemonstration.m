clear
clc
close all

% Generate a synthetic fHR sequence
x = zeros(1, 300);
x(50) = x(50) + 1000;
x(51) = x(51) - 1000;
x = cat(2, x, 15.0*ones(1, 50));
x = cat(2, x, zeros(1, 180));
x = cat(2, x, -28.0*ones(1, 30));
x = cat(2, x, zeros(1, 80));
x = cat(2, x, linspace(0, 35.0, 80));
x = cat(2, x, 35.0*ones(1, 70));
x = cat(2, x, linspace(35.0, 0.0, 50));
x = cat(2, x, zeros(1, 100));
x = cat(2, x, linspace(0, -37.0, 60));
x = cat(2, x, -37.0*ones(1, 150));
x = cat(2, x, linspace(-37.0, 0.0, 25));
x = cat(2, x, zeros(1, 500));
x = cat(2, x, 3.0*ones(1, 50));
x = cat(2, x, zeros(1, 100));
x = cat(2, x, linspace(0, 36.0, 500));
x = cat(2, x, 36.0*ones(1, 1400));
x = cat(2, x, linspace(36.0, 25.0, 250));
x = cat(2, x, 25.0*ones(1, 500));
x = cat(2, x, linspace(25.0, 10.0, 350));
x = cat(2, x, 10.0*ones(1, 1400));

x = x + 130; 

y = x + randn(size(x));
z = LPFilter(y, 0.03);
T = length(z);

% Find the baseline
baseline_est1 = median(z);
I = find(abs(z - baseline_est1) <= 15.0); % 15 BPM is the threshold for most accel/decel definitions
zz = z(I);
baseline_est2 = median(zz);

% Find the baseline
wlen = 1000;
BASELINE_RES = 5.0;
baseline_est3 = BaseLine1(z, wlen, 'md');

baseline_est4 = quant(baseline_est3, BASELINE_RES);

lbl = {};
figure
hold on
% plot(x);
plot(z, 'linewidth', 2); lbl = cat(2, lbl, {'fHR'}); 
plot(baseline_est1(ones(1, T)), 'linewidth', 2); lbl = cat(2, lbl, {'Global median'});
plot(baseline_est2(ones(1, T)), 'linewidth', 2); lbl = cat(2, lbl, {'Two-run global median'});
plot(baseline_est3, 'linewidth', 2); lbl = cat(2, lbl, {'Running median'});
plot(baseline_est4, 'linewidth', 2); lbl = cat(2, lbl, {['Running median quantized to ' num2str(BASELINE_RES) ' BPM']});
grid
legend(lbl, 'interpreter', 'none');
xlabel('Beat index');
ylabel('Beats per minute');