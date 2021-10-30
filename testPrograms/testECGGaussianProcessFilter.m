% An ECG denoiser based on a data driven Gaussian Process model of the ECG samples and noise

close all
clear
clc

% Load data
datafilepath = '../../../../DataFiles/PTBDataset/Physionet.org/files/ptbdb/1.0.0/';
directory_list = dir([datafilepath 'patient*']);
filelist = dir(fullfile([datafilepath, '**/*lr*.mat']));  % get list of all mat files
fs = 1000.0;

% Baseline wander removal filter
w1 = 0.72; % First stage baseline wander removal window size in seconds
w2 = 0.87; % Second stage baseline wander removal window size in seconds
BASELINE_REMOVAL_APPROACH = 'BP'; %'MDMN';
for k = 30 : 30%length(filelist)
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val;
    data = data(:, 1 : round(10*fs));
    
    switch(BASELINE_REMOVAL_APPROACH)
        case 'BP'
            data = data - LPFilter(data, 5.0/fs);
            data = LPFilter(data, 80.0/fs);
        case 'MDMN'
            wlen1 = round(w1 * fs);
            wlen2 = round(w2 * fs);
            for jj = 1 : size(data, 1)
                bl1 = BaseLine1(data(jj, :), wlen1, 'md');
                data(jj, :) = data(jj, :) - BaseLine1(bl1, wlen2, 'mn');
            end
    end
    
    ch = 2; % the desired channel
    sig = data(ch, :);
    SNR_pre_set = 15; % the desired input SNR
    sd = sqrt(var(sig) / 10^(SNR_pre_set/10));
    x = sig + sd * randn(size(sig));
    
    f0 = 1.2; % approximate heart rate (in Hz) used for R-peak detection
    peaks = PeakDetection(sig, f0/fs);                  % peak detection
    
    I = find(peaks);
    t = (0 : length(x) - 1)/fs;
    
    GPfilterparams.bins = 300; % number of phase domain bins
    GPfilterparams.NOISE_VAR_EST_METHOD = 'AVGLOWER';
    GPfilterparams.avg_bins = 3;
    GPfilterparams.SMOOTH_PHASE = 'GAUSSIAN';
    GPfilterparams.gaussianstd = 1.0;
    GPfilterparams.plotresults = 1;
    [data_posterior_est, data_prior_est] = ECGGaussianProcessFilter(x, peaks, GPfilterparams);   
    
    SNR_pre = 10 * log10(mean(sig.^2) / mean((x - sig).^2));
    SNR_prior = 10 * log10(mean(sig.^2) / mean((data_prior_est - sig).^2));
    SNR_posterior = 10 * log10(mean(sig.^2) / mean((data_posterior_est - sig).^2));
    
    figure;
    plot(t, x);
    hold on;
    plot(t(I), sig(I),'ro');
    grid
    xlabel('time (s)');
    ylabel('Amplitude');
    title('Noisy ECG and the detected R-peaks');
    set(gca, 'fontsize', 16)

    figure
    plot(t, x)
    hold on
    plot(t, data_posterior_est, 'linewidth', 2)
    plot(t, data_prior_est, 'linewidth', 2)
    plot(t, sig, 'linewidth', 2)
    grid
    legend('Noisy ECG', 'ECG posterior estimate', 'ECG prior estimate', 'Original ECG')    
    xlabel('time (s)');
    ylabel('Amplitude');
    title('Filtering results');
    set(gca, 'fontsize', 16)
    
    disp('Filtering performance:');
    disp([' Input SNR (desired) = ' num2str(SNR_pre_set), 'dB']);
    disp([' Input SNR (actual) = ' num2str(SNR_pre), 'dB']);
    disp([' Output SNR (prior-based) = ' num2str(SNR_prior), 'dB']);
    disp([' Output SNR (posterior-based) = ' num2str(SNR_posterior), 'dB']);
end

