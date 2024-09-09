% test IIR notch filter in C++ and in Matlab
clear;
close all;
x = load('testinput.txt')';
y = load('testoutput.txt')';

hparams = fopen('notchfilterparams.txt','r');
params = fscanf(hparams,'%g');

f0 = params(1);
fc = params(2);
fs = params(3);
fclose(hparams);

[b,a] = iirnotch(f0/(fs/2),fc/(fs/2));  
% s = filter(b,a,x);
% s = TrimmedFilter(x,length(x),'median',50)';
s = nonlinear_filter(x,'trmean',50,10,10,ones(50,1)/sum(ones(50,1)));
% s = iir_filter(b,a,x);

figure;
hold on;
plot(x(1:end-24));
plot(y(25:end),'r');
plot(s(25:end),'m');
% plot(s(1:end-24),'m');
grid;

figure;
hold on;
psd(x,1024,fs);
psd(y,1024,fs);
psd(s,1024,fs);

figure;
plot(abs(y(25:end)-s(1:end-24)));
% plot(20*log10(abs(y(25:end)-s(1:end-24))));
grid;