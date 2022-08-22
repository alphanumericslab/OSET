% Matlab script that uses the TikhonovRegularization.m;
% It reads a .csv file produced by testTikhonovRegularization.ipynb that
% contains a noisy ecg and it's corresponding reconstruction obtained by 
% the python function TikhonovRegularization.m . It calls the
% function TikhonovRegularization.m and adds as the third column the
% corresponding reconstructed ecg. The two are compared numerically in 
% testTikhonovRegularization.ipynb



datafilepath = fullfile('dataTikhonov.csv')
data = readtable(datafilepath);
ecg = data.ecg1;
%plot(ecg)

DiffOrderOrFilterCoefs = 4;
Lambda = 0.3;
tic
y = TikhonovRegularization(ecg', DiffOrderOrFilterCoefs, Lambda);
toc

% plot(ecg)
% hold on
% plot(y)
% size(y)

data.(size(data, 2)+1) = y';
data.Properties.VariableNames{numOfColumn+1} = 'ecg_hat_RS1';
%data
writetable(data, 'dataTikhonov.csv')
