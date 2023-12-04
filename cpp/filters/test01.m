% test IIR filter in C++ and in Matlab
clear;
close all;
x = load('testinput.txt')';
y = load('testoutput.txt')';

hcoefs = fopen('filtercoefs.txt','r');
data = fscanf(hcoefs,'%g');

M = data(1);
b = data(2:M+1);
N = data(M+2);
a = data(M+3:M+2+N);
fclose(hcoefs);


s = filter(b,a,x);

figure;
hold on;
% plot(x);
plot(y,'r');
plot(s,'mo');
grid;

figure;
plot(20*log10(abs(y-s)));
grid;