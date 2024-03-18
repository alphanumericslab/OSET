% Test script for removing EOG artifacts from EEG data using the
% iterative subspace denoising algorithm developed by Sameni and Gouy-Pailler (2014).
%
% Revision History:
%   2011: First release.
%   2020: Replaced obsolete psd() with pwelch() for power spectral density estimation.
%   2024: Renamed from deprecated version EOGRemoval().
%
% Reference:
%   Sameni, R. and Gouy-Pailler, C. (2014). An iterative subspace denoising
%   algorithm for removing electroencephalogram ocular artifacts. Journal of
%   neuroscience methods, 225, 97-105.
%
%   Sameni, R., Jutten, C., & Shamsollahi, M. B. (2010). A Deflation
%   Procedure for Subspace Decomposition. In IEEE Transactions on Signal
%   Processing (Vol. 58, Issue 4, pp. 2363–2374). Institute of Electrical
%   and Electronics Engineers (IEEE).
%   https://doi.org/10.1109/tsp.2009.2037353
%
% Reza Sameni, 2011-2024
% The Open-Source Electrophysiological Toolbox (OSET): https://github.com/alphanumericslab/OSET

clear
close all;
clc;

%% Select and load dataset
CASE_STUDY = 1; % select case study or add new ones
switch CASE_STUDY
    case 1
    load EEGdata
    fs = 250;
    data = data(fs*0+1:fs*20,1:25)';
    data = .1*(data - mean(data,2)*ones(1,size(data,2)))./(std(data,[],2)*ones(1,size(data,2))); % normalize data
    EEG = data([1:23 25],:);
    EOG = data(24,:);

    case 2
    load EEGdata2
    fs = 250;
    data = (data - mean(data,2)*ones(1,size(data,2)))./(std(data,[],2)*ones(1,size(data,2))); % normalize data
    EEG = data(1:end-1,:);
    EOG = data(end,:);
end

% Check noise-robustness (for research only)
% % % A = randn(50,size(EEG,1));
% % % EEG = A*EEG + .001*randn(50,size(EEG,2));

%% Set eog_canceller_nsca parameters
twlen = 0.3; % Time window length in seconds for EOG power envelope calculation
num_itr = 3; % Number of iterations for the artifact removal process
th_method = 'envelope-prctile'; % Threshold method for EOG segment detection. Options: 'absolute', 'envelope-mean-fraction', 'envelope-median-fraction', 'envelope-prctile'
th_value = 75.0; % Threshold value for the chosen method. Interpretation depends on th_method. In 'envelope-prctile' mode, this is the envelope amplitude percentile above which is EOG 
num_comp_est_method = 'opt'; % Method for estimating the number of components to denoise in each iteration. Options: 'fixed' or 'opt'
shrinkage_params = struct; % Initialize shrinkage parameters structure
switch num_comp_est_method % Configure parameters based on the component estimation method
    case 'fixed'
        % Use a fixed number of channels based on the largest generalized eigenvalues
        shrinkage_params.num_noise_gevds = 1;
    case 'opt'
        % Dynamically estimate the EOG dimensions embedded in background EEG noise
        % Additional parameters can be set here if required for dynamic estimation
end

shrinkage_method = 'wden'; % Choose the method for EOG component shrinkage. Options: 'remove', 'wden' for wavelet denoising
switch shrinkage_method % Configure shrinkage parameters based on the chosen method
    case 'remove'
        % Simply remove components without filtering, reducing dimensionality
        % This mode expects shrinkage_params.num_noise_gevds to be set
    case 'wden'
        % Wavelet denoising parameters
        shrinkage_params.TPTR = 'heursure'; % Threshold selection rule
        shrinkage_params.SORH = 's';        % Soft or hard thresholding
        shrinkage_params.SCAL = 'mln';      % Level-dependent threshold scaling
        shrinkage_params.NUM = 5;           % Wavelet decomposition level
        shrinkage_params.WNAME = 'sym5';    % Wavelet name
end
adapt_reference = 'top-few-gevd'; % Adaptation method for the reference EOG signal across iterations. Options: 'no', 'top-gevd', 'top-few-gevd'

flagplot = 3; % Plotting flag to control the verbosity of output plots. Try 1, 2 and 3

%% Call the EOG artifact removal function with configured parameters
EEG_cleaned = eog_canceller_nsca(EEG, EOG, fs, twlen, num_itr, th_method, th_value, num_comp_est_method, shrinkage_method, shrinkage_params, adapt_reference, flagplot);
% EEG_cleaned = EOGRemoval(EEG, EOG, round(twlen*fs), 1.0, 3, num_itr, 1, 0);

%% plot the results
t = (0:length(EEG)-1)/fs;
L = 4;
ch = 1:size(EEG,1);
for i = 1:length(ch)
    if(mod(i,L)==1)
        figure;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(t,EEG(ch(i),:),'k');
    hold on;
    plot(t,EEG_cleaned(ch(i),:),'color',.6*ones(1,3));
    % % %     grid;
    set(gca,'Box','On','FontSize',16);
    if (i<L)
        set(gca,'XTickLabel',[]);
    end
    ylabel(['EEG',num2str(i)]);
    axis tight
end
xlabel('time(s)','FontSize',16);

ch = 1; % Select the channel for spectral analysis
wlen = round(fs*1.0); % spectral window length
woverlap = round(fs*0.25); % spectral window overlap
nfft = 1000; % number of DFT points
figure
pwelch(EEG(ch,:), hamming(wlen), woverlap, nfft, fs);
hold on
pwelch(EEG_cleaned(ch,:), hamming(wlen), woverlap, nfft, fs);
legend('Original', 'After EOG Removal');
title(['EEG spectrum before and after EOG cancellation in channel ', num2str(ch)]);
set(gca, 'fontsize', 16)

EEG = (EEG - mean(EEG,2))./(std(EEG,[],2)*ones(1,size(EEG,2)));
n = 1:size(EEG,1);
egs = sort(log(eig(cov(EEG'))),'descend');
figure;
hold on;
plot(n,egs,'linewidth',3);
plot(n,egs,'ro','linewidth',3);
grid
set(gca,'fontsize',16,'box','on');
xlabel('n');
ylabel('$\lambda_n$ (dB)','interpreter','latex');
axis tight
title('Eigenvalue distribution of the EEG channels');
