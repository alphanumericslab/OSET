% test stability issues od the Kalman notch filter
% Reza Sameni, Copyright 2018

close all;
clear;
load('SampleECG1.mat'); %data = data';
% load('SampleECG2.mat'); data = data(1:15000,6)';
fs = 1000;
f0 = 50;

data = data - LPFilter(data,.7/fs);

n = (0:length(data)-1);
x = data + (.015 + .001*sin(2*pi*n*.1/fs)).*sin(2*pi*n*f0/fs)/std(data);



gamma = 1;
% [y1,y,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs,1.5e-8*max(abs(x)),var(data),gamma);
[y1,y,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs,1e-6*max(abs(x)),var(data),gamma);

% y = KFNotch(y,149.5,fs);
% y = KFNotch(y,249.2,fs);
% y = KFNotch(y,349,fs);
% y = KFNotch(y,448,fs);

Pbar = squeeze(Pbar(1,1,:));
Phat = squeeze(Phat(1,1,:));
PSmoothed = squeeze(PSmoothed(1,1,:));

I = 1:10*fs;

t = (0:length(x)-1)/fs;

% % % figure;
% % % plot(t,data);
% % % grid;
% % % title('original');

figure;
hold on;
plot(t(I),x(I));
plot(t(I),y(I)-.8,'r');
plot(t(I),data(I)-1.6,'g');
grid;
legend('Noisy','PL removed','Original');
xlabel('time(s)','FontSize',16);
ylabel('Amplitude (mV)','FontSize',16);
set(gca,'Box','On','FontSize',16);
set(gcf,'Position',[438 350 860 420]);


figure;
plot(t,Kgain');
grid;
title('KF gains');

% % % figure;
% % % hold on;
% % % plot(Pbar,'b');
% % % plot(Phat,'r');
% % % plot(PSmoothed,'g');
% % % grid;

figure;
hold on;
psd(x,10*fs,fs);
psd(y,10*fs,fs);
psd(data,10*fs,fs);
legend('Noisy','Denoised','Original');


c = cos(2*pi*f0/fs);
k = Kgain(:,end);
a = 1/(1-k(1));
figure;
freqz([1 -2*c 1],[a -a*k(2)-2*c 1],10*fs,fs);


clear H;
c = cos(2*pi*f0/fs);
f = (49.99:.00001:50.01);
for i = 1:5:1000
    k = Kgain(:,i);
    a = 1/(1-k(1));
    H(i,:) = freqz([1 -2*c 1],[a -a*k(2)-2*c 1],f,fs);
end

figure;
hold on;
plot(f,20*log10(abs(H))');
grid;



