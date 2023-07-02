clear
close all
clc

% Fixed-Order variable Lambda
DiffOrder = 2; % smoothness constraint order greater than 1
% lambda = [1e-2 : 1e-2 : 100 500 1000 10000]; % smoothness factor
lambda = [1e-4 2.5e-4 5e-4 1e-3 2.5e-3 5e-3 7.5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1:1e-1:100 500 5000 1e4 2e4 1e5]; % smoothness factor
% lambda = logspace(-3, 10, 1000);
h = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder);
phi = conv(h,h);
N = 512;
r = zeros(length(lambda), 2*DiffOrder);
H = zeros(N, length(lambda));
W = zeros(N, length(lambda));
for i = 1 : length(lambda),
    gg = phi*lambda(i);
    midindex = (length(phi)+1)/2;
    gg(midindex) = gg(midindex) + 1;
    r(i,:) = roots(gg);
    [H(:, i), W(:, i)] = freqz(sum(gg),gg, N);
end

% Fixed-Lambda variable order
DiffOrder = [2 4 6 8 10 12 14]; % smoothness constraint order greater than 1
llambda = 1; % smoothness factor

N = 1024;
HH = zeros(N, length(DiffOrder));
WW = zeros(N, length(DiffOrder));
for i = 1 : length(DiffOrder),
    hh = diff([zeros(1, DiffOrder(i)) 1 zeros(1, DiffOrder(i))], DiffOrder(i));
    phiphi = conv(hh,hh);
    gg = phiphi*llambda;
    %     gg = phiphi*llambda^(DiffOrder(i));
    midindex = (length(phiphi)+1)/2;
    gg(midindex) = gg(midindex) + 1;
    [HH(:, i), WW(:, i)] = freqz(sum(gg),gg, N);
end

theta = -1.1*pi:.1:pi;

figure
hold on
% plot(r,'bx', 'linewidth', 2);
plot(r,'b', 'linewidth', 2);
plot(cos(theta),sin(theta),'r--', 'linewidth', 3);
axis([-1 1.8 -1.2 1.2]);
axis equal
set(gca,'box','on');
set(gca,'FontSize',18);
grid on;
text(-.45, .1, '\lambda_0 = 0', 'FontSize', 18);
text(1.1, 0, '\lambda_0 = \infty', 'FontSize', 18);
text(.25, .38, '\rightarrow', 'FontSize', 18);
text(.25, -.38, '\rightarrow', 'FontSize', 18);
text(1.55, 0.9, '\downarrow', 'FontSize', 18);
text(1.55, -1, '\uparrow', 'FontSize', 18);
% % % text(1., 1.05, '\lambda_0 = 0', 'FontSize', 18);
% % % text(1., -1.05, '\lambda_0 = 0', 'FontSize', 18);

figure
fg2 = plot(W(:, 1)/pi, (abs(H(:,[1, 2, 3, 4, 5, 6, 8, 12, 20, 30, 40, 60, 100, 200, 300, 500, 1002, 1003, 1004, end]))), 'linewidth', 2);
grid
set(gca,'FontSize',18);
axis([0 1 -.05 1.05]);
xlabel('Normalized frequency (\omega/\pi)');
ylabel('Amplitude');
set(gca,'box','on');

figure
fg3 = plot(WW(:, 1)/pi, abs(HH), 'linewidth', 2);
grid
set(gca,'FontSize',18);
axis([0 1 -.05 1.05]);
xlabel('Normalized frequency (\omega/\pi)');
ylabel('Amplitude');
set(gca,'box','on');

% figuresize(3 , 3);

% figure
% plot(W(:, 1)/pi, (angle(H(:,[1, 2, 3, 4, 5, 6, 8, 12, 20, 30, 40, 60, 100, 200, 300, 500, 1002, 1003, 1004, end]))), 'linewidth', 2);
% grid
% set(gca,'FontSize',18);
% % axis([0 1 0 1]);
% xlabel('Normalized frequency (2f/f_s)');
% ylabel('Phase (rad)');
% set(gca,'box','on');

% % % pos = get(gca,'Position');
% % % arrowObj = annotation('arrow',[.01 .012] + .3,[.228 .32] + pos(2));

% h = quiver(0.01,0.15,0.1*cos(pi/3.5),0.1*sin(pi/3.5));
% set(h, 'color', 'k', 'linewidth', 1.5, 'MarkerSize', 12, 'MaxHeadSize', 5);
% set(gca, 'XLim', [1 10], 'YLim', [1 10]);

% set(arrowObj, 'Units', 'centimeters');
% set(arrowObj, 'Position', [1 1 3 5]);
% set(gca,'LineWidth',2);

% % % num = 1;
% % % den = phi;
% % %
% % % sys = tf(num, den);
% % % zero = roots(num);
% % % pole = roots(den);
% % %
% % % figure
% % % rlocus(-sys, lambda);
% % %
% % % [R K] = rlocus(-sys, lambda);
% % % figure;
% % % hold on;
% % % for i = 1:size(R,1)
% % %     plot(real(R(i,:)),imag(R(i,:)),'linewidth',3);
% % % end
% % % for i = 1:length(pole)
% % %     plot(real(pole(i)),imag(pole(i)),'kx','linewidth',2,'markersize',12);
% % % end
% % % for i = 1:length(zero)
% % %     plot(real(zero(i)),imag(zero(i)),'ko','linewidth',2,'markersize',12);
% % % end
% % % xlabel('Real','fontsize',14);
% % % ylabel('Imaginary','fontsize',14);
% % % set(gca,'fontsize',14);
% % % set(gca,'box','on');

%
% gg = phi*SmoothnessFactor;
% midindex = (length(phi)+1)/2;
% gg(midindex) = gg(midindex) + 1;
% r = roots(gg);
% r_abs = abs(r);
% I_causal = r_abs < 0.9999;
% r_causal = r(I_causal);
% p_causal = poly(r_causal);
% y = filtfilt(sum(p_causal), p_causal, x);
%
% close all;
% clear;
% fs = 256;
% f0 = 50;
%
