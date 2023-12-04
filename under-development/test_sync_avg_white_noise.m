% clone thid first: https://github.com/alphanumericslab/OSET/tree/staging

clear
close all
clc

fs = 250;
x = lp_filter_zero_phase(randn(1, 100000), 1.5/fs);
x = x - lp_filter_zero_phase(x, 1.49/fs);

sign_changes = x(1:end-1) < 0 & x(2:end) >= 0; % case I
% sign_changes = x(1:end-1) .* x(2:end) <= 0; % case II
I_sign_changes = find(sign_changes);

event_width = round(5.0*fs);
if mod(event_width, 2) == 0
    event_width = event_width + 1;
end
y = event_stacker(x, I_sign_changes, event_width);

[mn, vr_mn, md, vr_md] = robust_weighted_average(y);

figure
plot(x)

figure
hold on
plot(y')
plot(mn, 'k', 'linewidth', 3)
plot(cumsum(mn), 'r', 'linewidth', 3)
grid
