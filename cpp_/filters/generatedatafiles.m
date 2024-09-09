clear;
close all;

%//////////////////////////////////////////////////////////////////////////
% generate input data file
x = cumsum(randn(10000,1));
x = x + sin(2*pi*(1:length(x))'*60/1000);

hin = fopen('testinput.txt','w');
fprintf(hin,'%14.10f\n',x);
fclose(hin);

%//////////////////////////////////////////////////////////////////////////
% generate IIR filter coefficients
b = hamming(150);
a = 1;

hcoefs = fopen('filtercoefs.txt','w');
fprintf(hcoefs,'%d\n',length(b));
fprintf(hcoefs,'%14.10f',b);
fprintf(hcoefs,'\n%d\n',length(a));
fprintf(hcoefs,'%14.10f',a);
fclose(hcoefs);

%//////////////////////////////////////////////////////////////////////////
% set notch filter parameters
f0 = 60;
fc = 1;
fs = 1000;

hcoefs = fopen('notchfilterparams.txt','w');
fprintf(hcoefs,'%14.10f \n %14.10f \n %14.10f \n',f0, fc, fs);
fclose(hcoefs);
