% function RunningXcorr for ECG phase wrapping analysis
% Reza Sameni, 2015
%
close all;
clear;
clc

data = load('ECG1.txt');
fs = 1000;

% data = load('16265m.mat'); data = data.val(1,:);
% fs = 128;

data = LPFilter(data - LPFilter(data,5/fs),50/fs);

% x = data(1001:11000);
% y = data(1001:11000);
x = data(1:5000);
y = data(601:1500);

x = x(:)';
y = y(:)';

r = x'*y;

% % % % x = randn(1,1000);
% % % % y = x;
% % % 
% % % m = 1:1:round(fs/.5);
% % % len = 10;
% % % % % % len2 = 5;
% % % 
% % % 
% % % if(length(x)>length(y))
% % %     y = [y zeros(1,length(x)-length(y))];
% % % elseif(length(y)>length(x))
% % %     x = [x zeros(1,length(y)-length(x))];
% % % end
% % % 
n = 1:length(x);
% % % r = zeros(length(m),length(x));
% % % for k = 1:length(m)
% % %     xx = x(1:end-m(k)+1);
% % %     yy = y(m(k):end);
% % %     
% % %     zz = xx.*yy;
% % %     
% % %     %         rr = filter(ones(1,len),len,zz);
% % %     %     rr = filter(hamming(len),1,zz);
% % %     rr = filtfilt(hamming(len),1,zz);
% % %     
% % %     % % %     e = filter(ones(1,round(len2)),round(len2),zz);
% % %     % % %     ee = sqrt(filter(ones(1,round(len2)),round(len2),zz.^2));
% % %     % % % %     rr = filter(ones(1,len),len,zz./ee);
% % %     % % %     rr = filter(ones(1,len),len,zz./sqrt(ee.^2-e.^2));
% % %     
% % %     r(k,:) = [rr zeros(1,m(k)-1)];
% % % end
% % % 
% % % % peaks = PeakDetection(x,1.5/fs);
% % % % I = find(peaks);
% % % % d = diff(I);
figure;
subplot(211);
plot(n/fs,x);
grid
axis tight
% subplot(312);
% plot(d);
% grid
subplot(212);
% image(r')
% r = r - min(min(r));
% r(r<0) = 0;
% r = log(abs(r)+1);
r = r - min(r(:));
r = r./max(r(:));
% image(n/fs,m/fs,255*r);
image(255*r');

% image(n/fs,m/fs,255*r/max(max(abs(r))));
hold on
plot(n/fs,4*x+.5,'r');

% figure
% mesh(n/fs,m,log(abs(r)));
% grid
% colormap cool

