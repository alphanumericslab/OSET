close all
clear
clc

% Convert all dat files in bash using this command: find ./*/ -type f -execdir wfdb2mat -r {} \;

datafilepath = '../../../../DataFiles/PTBDataset/Physionet.org/files/ptbdb/1.0.0/';
directory_list = dir([datafilepath 'patient*']);

filelist = dir(fullfile([datafilepath, '**/*lr*.mat']));  % get list of all mat files

fs = 1000.0;
w1 = 0.72; % First stage baseline wander removal window size in seconds
w2 = 0.87; % Second stage baseline wander removal window size in seconds
BASELINE_REMOVAL_APPROACH = 'BP'; %'MDMN';
for k = 30 : 30%length(filelist)
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val;
    
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
    
    ch = 1;
    sd = 0.5 * std(data(ch, :))
    x = data(ch, :) + sd * randn(size(data(ch, :)));

    f0 = 1.2;
    peaks = PeakDetection(data(ch, :), f0/fs);                  % peak detection
    
    I = find(peaks);
    t = (0 : length(x) - 1)/fs;
    figure;
    plot(t, x);
    hold on;
    plot(t(I), data(ch, I),'ro');
    grid
    
    [phase, phasepos] = PhaseCalculation(peaks);     % phase calculation
    
    teta = 0;                                       % desired phase shift
    pphase = PhaseShifting(phase,teta);             % phase shifting
    
    bins = 250;                                     % number of phase bins
    [ECG_mean, ECG_std, meanphase] = MeanECGExtraction(x, pphase, bins, 1); % mean ECG extraction
    
    % Method 1: min std
    %     noise_std_est = min(ECG_std)
    
    % Method 2: average of the smallest std
    avg_bins = 5;
    ECG_std_up_sorted = sort(ECG_std);
    %     noise_std_est = median(ECG_std_up_sorted(1 : avg_bins))
    noise_std_est = sqrt(mean((ECG_std_up_sorted(1 : avg_bins).^2)))
    
    % Method 3: percentiles
    %     p = 0.5; % the desired percentile
    %     noise_std_est = prctile(ECG_std, p)
    
    ECG_intrinsic_std = sqrt(max(0, ECG_std.^2 - noise_std_est^2));
    
    lambda = 1.0;
    ECG_intrinsic_std_smoothed = TikhonovRegularization(ECG_intrinsic_std, 2, lambda);
    
    figure
    errorbar(meanphase, ECG_mean, ECG_std);
    grid
    
    figure
    plot(meanphase, ECG_std);
    hold on
    plot(meanphase, ECG_intrinsic_std);
    plot(meanphase, ECG_intrinsic_std_smoothed);
    grid
    legend('ECG + Noise STD', 'ECG STD', 'ECG Smoothed STD');
    
    figure
    plot(sort(ECG_std));
    grid
end

