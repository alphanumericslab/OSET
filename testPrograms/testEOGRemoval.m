% test file for removing EOG artifacts from the EEG
% Reza Sameni, Copyright 2011
% Created 2011
% Modified June 2018
%
% Ref:
% Sameni, R. and Gouy-Pailler, C. (2014). An iterative subspace denoising
% algorithm for removing electroencephalogram ocular artifacts. Journal of 
% neuroscience methods, 225, 97-105.


clear
close all;
clc;

load EEGdata
fs = 250;
data = data(fs*0+1:fs*20,1:25)';
data = .1*(data - mean(data,2)*ones(1,size(data,2)))./(std(data,[],2)*ones(1,size(data,2))); % normalize data
EEG = data([1:23 25],:);
EOG = data(24,:);

% % % load EEGdata2
% % % fs = 250;
% % % data = (data - mean(data,2)*ones(1,size(data,2)))./(std(data,[],2)*ones(1,size(data,2))); % normalize data
% % % EEG = data(1:end-1,:);
% % % EOG = data(end,:);

% % % A = randn(50,size(EEG,1));
% % % EEG = A*EEG + .001*randn(50,size(EEG,2));

twlen = 0.3; % in seconds
wlen = round(twlen*fs);
th = 1;
M = 3;
L = 8;

y = EOGRemoval(EEG,EOG,wlen,th,M,L,1,0);

% % % I = fs*14+1:fs*29;
% % % EEG = EEG(:,I);
% % % y = y(:,I);

t = (0:length(EEG)-1)/fs;
L = 4;
% % % % ch = [1 5 15 23];
ch = 1:size(EEG,1);
for i = 1:length(ch),
    if(mod(i,L)==1)
        figure;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(t,EEG(ch(i),:),'k');
    hold on;
    plot(t,y(ch(i),:),'color',.6*ones(1,3));
    % % %     grid;
    set(gca,'Box','On','FontSize',16);
    if (i<L)
        set(gca,'XTickLabel',[]);
    end
    ylabel(['EEG',num2str(i)]);
    axis tight
end
xlabel('time(s)','FontSize',16);

figure
% h = spectrum.welch;
% psd(h,EEG(ch(i),:),'Fs', fs, 'NFFT', 1024);
psd(EEG(ch(1),:), 1000, fs);
hold on
% psd(h,y(ch(i),:),'Fs', fs, 'NFFT', 1024);
psd(y(ch(1),:), 1000, fs);
legend('Original', 'After EOG Removal');


% % % EEG = (EEG - mean(EEG,2)*ones(1,size(EEG,2)))./(std(EEG,[],2)*ones(1,size(EEG,2)));

% % % n = 1:size(EEG,1);
% % % egs = sort(log(eig(cov(EEG'))),'descend');
% % % figure;
% % % hold on;
% % % plot(n,egs,'linewidth',3);
% % % plot(n,egs,'ro','linewidth',3);
% % % grid
% % % set(gca,'fontsize',16,'box','on');
% % % xlabel('n');
% % % ylabel('$\lambda_n$ (dB)','interpreter','latex');
% % % axis tight
