% Test the instantaneous correlation of different signals
% Reza Sameni

clc
clear;
close all;

N = 1000;
fs = 1000;
n = 0:N-1;
x = sin(2*pi*n*10/fs);
y = sin(2*pi*n*33/fs);
z = [x y];
zz = z + .05*randn(size(z));
zz = [zz randn(1,1000) zz ];

% % % load sampledata x fs
% % % zz = x(1,:);

tau = 0:199;
t = (0:length(zz)-1)/fs;
r = InstCorr(zz,tau);
mesh(t/fs,tau,r');xlabel('t');ylabel('tau');

% % % figure
% % % plot((0:2*N-1)/fs,zz);
% % % grid;
% % % figure
% % % mesh(t/fs,tau,r');xlabel('t');ylabel('tau');

% % % plot(sum(r,1))
% % % plot(sum(r,2));grid
% % % [r t tau] = instcorr(zz(end/2-600:end/2+600));
% % % plot(sum(r,2));grid
% % % ww = [zz randn(1,1000) zz];
% % % [r t tau] = instcorr(ww(end/2-800:end/2+800));
% % % plot(sum(r,2));grid

figure;
subplot(211)
plot(t, zz);grid
subplot(212)
plot(t, sum(r,2));grid
xlabel('time(s)');
% % % [r t tau] = instcorr(ww(end/2-1000:end/2+1000));
% % % figure
% % % plot(sum(r,2));grid
% % % [r t tau] = instcorr(ww(end/2-2000:end/2+2000));
% % % plot(sum(r,2));grid
% % % [r t tau] = instcorr(ww(end/2-2000:end/2+3000));
% % % [r t tau] = instcorr(ww(end/2-3000:end/2+3000));
% % % plot(sum(r,2));grid
% % % clear;
% % % close all;
% % % N = 5000;
% % % fs = 1000;
% % % n = 0:N-1;
% % % x = sin(2*pi*n*10/fs);
% % % y = sin(2*pi*n*33/fs);
% % % z = [x y];
% % % zz = z + .05*randn(size(z));
% % % plot((0:2*N-1)/fs,zz)
% % % [r t tau] = instcorr(zz);
% % % mesh(t/fs,tau,r');xlabel('t');ylabel('tau');