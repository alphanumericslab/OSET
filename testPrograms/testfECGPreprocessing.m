% test fetal ECG preprocessing
% Reza Sameni
% July 2018

clc
clear;
close all;

load('FOETAL_ECG.dat'); data = FOETAL_ECG(:,2:end)'; time = FOETAL_ECG(:,1)'; clear FOETAL_ECG; fs = 250;

plot(data');
PlotECG(data, 4, 'b', fs)