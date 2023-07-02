clear
close all
clc

% load ECG1
% s = data;

load ECG2
% s = data([1 4 5 6 12 15], round(0.3*fs):round(6.3*fs));
% s = data([1 3 4 5 6 7 8 12 15], round(0.3*fs):round(6.3*fs));
s = data(:, round(0.3*fs):round(6.3*fs));
s = LPFilter(s - LPFilter(s, 1.0/fs), 100.0/fs);

DiffOrder = 2; % smoothness constraint order greater than 1
snr = 15; % dB

s = s - mean(s, 2)*ones(1, size(s,2));
s = s./(std(s, [], 2)*ones(1, size(s,2)));

nvar = var(s, [], 2)/10^(snr/10);
noise = sqrt(nvar(:, ones(1, size(s, 2)))).*randn(size(s));
x = s + noise;

% h = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder);
load DiffCoefs
L = length(h);
N = size(x, 2);
% D = sparse(toeplitz([h(1) zeros(1, N - L)], [h zeros(1, N - L)])');
D = sparse(toeplitz([h(1) zeros(1, N - L + 2)], [h zeros(1, N - L)])');

a = (x*D);
A = a*a';
B = x*x';

% % lambda = 100000;
% % a = x*((lambda*(D*D') + eye(N))\D);
% x = x'; % make the dimensions TxN
% a = lambda*x'*((D*D' + lambda*eye(N))\D);
% A = a*a';
% % % x = x'; % make the data back to NxT

[V,d] = eig(A,B);
% [V,d] = eig(B);
d = diag(d);
y = V'*x;

% gamma = .1;
% N = size(D, 2);
% X = x';
% aa = X'*((gamma*(D'*D) + eye(N))\D');
% AA = aa*aa';
% BB = X'*X;
%
% [VV,dd] = eig(AA,BB);
% diag(dd)
% Y = X*VV;
% yy = Y';

yy = y;
yy(4 : end, :) = 0;
xx = pinv(V')*yy;

% PlotECG(x, 6, 'b', fs);
% PlotECG(y, 6, 'r', fs);
% PlotECG(yy, 5, 'm', fs);


% ch = 1;
% figure
% hold on
% plot(y(ch,:));
% plot(yy(ch,:),'r');
% plot(x(ch,:)*D,'k');
% grid

% % % x_ = (x - mean(x,2)*ones(1,size(x,2)))./(std(x,[],2)*ones(1,size(x,2))); % normalize data
% % % [lambda AIC MDL NEW ENSTh ENS fhand] = EstimateDimension(x_, var(x_(1, :)), 1);
[lambda AIC MDL NEW ENSTh ENS fhand] = EstimateDimension(y, nvar(1), 1);

% figure(fhand(3));
% hold on
% dd = sort(d,'ascend');
% plot(log(dd))

bias = .19;
L = size(x, 1);
t = (0 : size(x, 2) - 1)/fs;
figure('Position', [166 105 680 707]);
hold on
for i = 1:L,
    scale = norm(s(i,:));
    plot(t,x(i,:)/scale - bias*i,'linewidth', 2, 'color', 0.0*ones(1, 3));
    %     plot(t,s(i,:)/scale - bias*i,'linewidth', 2, 'color', 0.0*ones(1, 3));
end
grid;
set(gca, 'YTickLabel', [])
set(gca, 'YTick', []);
set(gca, 'fontsize', 16);
axis tight
xlabel('time(s)', 'fontsize', 16);
set(gca, 'Box', 'on');

figure('Position', [166 105 680 707]);
hold on
for i = 1:L,
    scale = norm(y(i,:));
    plot(t,y(i,:)/scale - bias*i,'linewidth', 2, 'color', 0.0*ones(1, 3));
end
grid;
set(gca, 'YTickLabel', [])
set(gca, 'YTick', []);
set(gca, 'fontsize', 16);
axis tight
xlabel('time(s)', 'fontsize', 16);
set(gca, 'Box', 'on');

figure('Position', [166 105 680 707]);
hold on
for i = 1:L,
    scale = norm(s(i,:));
    %     plot(t,x(i,:)/scale - bias*i,'linewidth', 2, 'color', 0.5*ones(1, 3));
    plot(t,xx(i,:)/scale - bias*i,'linewidth', 2, 'color', 0.0*ones(1, 3));
end
grid;
set(gca, 'YTickLabel', [])
set(gca, 'YTick', []);
set(gca, 'fontsize', 16);
axis tight
xlabel('time(s)', 'fontsize', 16);
set(gca, 'Box', 'on');


nn = 1:L;
dth = d(4);
figure
hold on
% stem(d/max(d), 'k', 'linewidth', 2);
% plot(nn, dth(ones(1,length(nn))), '--','linewidth',4,'color',.5*ones(1,3));
plot(nn, 5*ones(1,length(nn)), '--','linewidth',4,'color',.5*ones(1,3));
plot(nn, d, 'k', 'linewidth', 3);
plot(nn, d,'.','markersize',30,'color',.2*ones(1,3)); % all channels should be normalized to have the same variance
grid
% set(gca,'yscale','log');
set(gca, 'fontsize', 16);
xlabel('$n$', 'fontsize', 16, 'interpreter', 'latex');
% ylabel('$\frac{l(n)}{l_{max}}$', 'fontsize', 16, 'interpreter', 'latex');
ylabel('$l_n$', 'fontsize', 16, 'interpreter', 'latex');
axis tight;
set(gca, 'Box', 'on');
% a = axis;
% a(1) = 1;
% a(2) = L;
% axis(a);
