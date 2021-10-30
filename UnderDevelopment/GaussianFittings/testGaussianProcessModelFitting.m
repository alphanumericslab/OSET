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
    
    ch = 2;
    sig = data(ch, :);
    sd = 1 * std(sig)
    x = sig + sd * randn(size(sig));
    
    f0 = 1.2;
    peaks = PeakDetection(sig, f0/fs);                  % peak detection
    
    I = find(peaks);
    t = (0 : length(x) - 1)/fs;
    figure;
    plot(t, x);
    hold on;
    plot(t(I), sig(I),'ro');
    grid
    
    [phase, phasepos] = PhaseCalculation(peaks);     % phase calculation
    
    teta = 0;                                       % desired phase shift
    pphase = PhaseShifting(phase,teta);             % phase shifting
    
    bins = 250;%median(diff(I));                                     % number of phase bins
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
    
    % Convert ECG phase into a (phase x time) matrix
    M = ECGPhaseToMatrix(pphase, bins);
    
    % smooth the matrix in time and phase
    wlen_time = 1;
    wlen_phase = 1;
    %     M_smoothed = filter2(ones(wlen_phase, wlen_time)/(wlen_phase * wlen_time), M);
    M_smoothed = imgaussfilt(M, wlen_time);
    
    % reconstruct the ECG
    s_mean = ECG_mean * M_smoothed;
    %s_var = (ECG_intrinsic_std.^2) * M_smoothed;
    s_var = (ECG_intrinsic_std_smoothed.^2) * M_smoothed;
    n_var = noise_std_est^2;
    n_mean = 0;
    
    x_reconstructed = (s_var .* (x - n_mean) + n_var * s_mean) ./ (s_var + n_var);
    
    
    SNR_pre = 10 * log10(mean(sig.^2) / mean((x - sig).^2))
    SNR_smean = 10 * log10(mean(sig.^2) / mean((s_mean - sig).^2))
    SNR_post = 10 * log10(mean(sig.^2) / mean((x_reconstructed - sig).^2))
    
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
    
    figure
    imagesc(M_smoothed)
    
    figure
    plot(x)
    hold on
    plot(sig, 'linewidth', 2)
    plot(x_reconstructed, 'linewidth', 2)
    plot(s_mean, 'linewidth', 2)
    grid
    legend('Noisy ECG', 'Original ECG', 'Reconstruced from phase domain mean ECG', 'smean')
    
    figure
    imagesc(M_smoothed' * diag(ECG_intrinsic_std) * M_smoothed);
end

