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
writetable(data, '/Users/mircea/Work/GitHub/OSET/Python3/testPrograms/dataTikhonov.csv')
