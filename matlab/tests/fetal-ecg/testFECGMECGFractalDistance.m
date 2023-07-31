% test fetal ECG fractal dimension calculation
% Reza Sameni, Copyright March 2006
% Modified June 2018
%

clear;
close all;
clc;
%load patient165_s0323lre; data = data(1:10000,4)';
load FOETAL_ECG.dat; data = FOETAL_ECG(:,4)';clear FOETAL_ECG;
figure;
plot(data');grid;
%plot(data(6:8,:)');grid;

N = length(data);
for delta = 1:400;
    cntr = 0;
    s = 0;
    for m = 1:N-delta,
        %n = [m-delta:-delta:1 m+delta:delta:N];
        n = [m+delta:delta:N];
        s = s + sum(sum(abs(data(:,m)*ones(1,length(n))-data(:,n)).^1));
        %s = s + sum(sum( abs(data(:,m)*ones(1,length(n))-data(:,n))./ ( abs((data(:,m)*ones(1,length(n))).*data(:,n))) ));
        cntr = cntr + length(n);
    end
    C(delta) = s/cntr;
end

figure;
plot(C);
grid;

% figure;
% psd(data,length(data),1000);
