clear
close all
clc

% load ECG1
% s = data;

load ECG2
s = data(:, 1:5*fs);%end);
s = LPFilter(s - LPFilter(s, 1.0/fs), 100.0/fs);

DiffOrder = 2; % smoothness constraint order greater than 1
snr = 5; % dB

nvar = var(s, [], 2)/10^(snr/10);
x = s + sqrt(nvar(:, ones(1, size(s, 2)))).*randn(size(s));

h = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder);
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
diag(d)
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
yy(5:end, :) = 0;
xx = pinv(V')*yy;

PlotECG(x, 5, 'b', fs);
PlotECG(y, 5, 'r', fs);
PlotECG(xx, 5, 'm', fs);

ch = 1;
figure
hold on
plot(y(ch,:));
% plot(yy(ch,:),'r');
% plot(x(ch,:)*D,'k');
grid