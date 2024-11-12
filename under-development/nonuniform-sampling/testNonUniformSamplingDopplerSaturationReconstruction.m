close all
clear
clc

% Convert all dat files in bash using this command: find ./*/ -type f -execdir wfdb2mat -r {} \;

datafilepath = '../../../../DataFiles/physionet.org/files/ptbdb/1.0.0/';
directory_list = dir([datafilepath 'patient*']);

filelist = dir(fullfile([datafilepath, '**/*lr*.mat']));  % get list of all mat files

w1 = 0.72; % First stage baseline wander removal window size in seconds
w2 = 0.87; % Second stage baseline wander removal window size in seconds
BASELINE_REMOVAL_APPROACH = 'BYPASS';%'BP'; %'MDMN';

datafilename = '/Users/rsameni/Documents/DataFiles/Acoustic_1D_Doppler_OneDrive/ultrasound_p55v1_2_725985530069.wav';
% data = load(datafilename);
% data = data.val;

[data0, fs] = audioread(datafilename, 'native');
data0 = double(data0(:)');

switch(BASELINE_REMOVAL_APPROACH)
    case 'BP'
        data = data0 - lp_filter_zero_phase(data0, 1.0/fs);
        data = lp_filter_zero_phase(data, 1000.0/fs);
    case 'MDMN'
        wlen1 = round(w1 * fs);
        wlen2 = round(w2 * fs);
        data = zeros(size(data0));
        for jj = 1 : size(data, 1)
            bl1 = baseline_sliding_window(data0(jj, :), wlen1, 'md');
            data(jj, :) = data(jj, :) - baseline_sliding_window(bl1, wlen2, 'mn');
        end
    otherwise
        data = data0;
end

% preserve_percentage = 0.6;
ITR = 10;
T = size(data, 2);
% N_keep = round(T * preserve_percentage);
% kept_indexes = randi(T, [1, N_keep]);

% fc = 1000.0;
% fc = linspace(750, 1000, ITR);
fc = 1000.0 + zeros(1, ITR);

th = 28000;%0.86;
kept_indexes = find(abs(data) < th);

ecg_org = data;%lp_filter_zero_phase(data, fc(1) / fs);
ecg_dc = mean(ecg_org);
ecg_zero_mean = ecg_org - ecg_dc;

ecg_rec = interp1(kept_indexes, ecg_zero_mean(kept_indexes), 1:T);%zeros(1, T);
ecg_rec(kept_indexes) = ecg_zero_mean(kept_indexes);

% load lpf_1khz_at_48khz;
% lpf_1khz_at_48khz = lpf_1khz_at_48khz / sum(lpf_1khz_at_48khz);

ecg_all = zeros(ITR, T);
for m = 1 : ITR
    ecg_all(m, :) = ecg_rec;
    ecg_rec(kept_indexes) = ecg_zero_mean(kept_indexes);

    ecg_rec_lp = lp_filter_zero_phase(ecg_rec, fc(m) / fs);
    % ecg_rec_lp = lp_filter_zero_phase(ecg_rec, fc / fs);
    % ecg_rec_lp = filtfilt(lpf_1khz_at_48khz, 1, ecg_rec);
    
    ecg_rec = ecg_rec_lp/sum(ecg_rec_lp.^2)*sum(ecg_rec.^2);
    % ecg_rec = ecg_rec_lp;%/sum(ecg_rec_lp.^2)*sum(ecg_rec.^2);
    %         hr_continuous = BPFilter5(hr_continuous, 10.0/fs, 35.0/fs, 3); % BPFilter5(x,fc,bw,order)
    disp(['iteration = ', num2str(m)])
end
ecg_all = ecg_all + ecg_dc;
ecg_rec = ecg_rec + ecg_dc;

lgnd = {};
t = (0 : T - 1) / fs;
figure
hold on
plot(t, data0, 'linewidth', 2); lgnd = cat(2, lgnd, {'original'});
% plot(t, ecg_org, '-*', 'markersize', 14, 'linewidth', 2); lgnd = cat(2, lgnd, {'original'});
%     plot(t, ecg_all');
plot(t, ecg_all(end, :), 'r', 'linewidth', 1); lgnd = cat(2, lgnd, {'reconstructed'});
%     plot(t, ecg_all(1, :), 'b', 'linewidth', 2); lgnd = cat(2, lgnd, {'first iteration'});
%     plot(t, (ecg_all(1, :) + ecg_all(end, :)) / 2.0, 'r', 'linewidth', 2);
% plot(t(kept_indexes), ecg_org(kept_indexes), 'ro', 'markersize', 16); lgnd = cat(2, lgnd, {'kept points'});
grid
legend(lgnd);
xlabel('time(s)');
ylabel('Amplitude');
set(gca, 'fontsize', 16)


