close all
clear
clc

% Convert all dat files in bash using this command: find ./*/ -type f -execdir wfdb2mat -r {} \;

datafilepath = '../../../../DataFiles/physionet.org/files/ptbdb/1.0.0/';
directory_list = dir([datafilepath 'patient*']);

filelist = dir(fullfile([datafilepath, '**/*lr*.mat']));  % get list of all mat files

fs = 1000.0;
w1 = 0.72; % First stage baseline wander removal window size in seconds
w2 = 0.87; % Second stage baseline wander removal window size in seconds
BASELINE_REMOVAL_APPROACH = 'BP'; %'MDMN';
for k = 1 : 5%length(filelist)
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val;

    switch(BASELINE_REMOVAL_APPROACH)
        case 'BP'
            data = data - LPFilter(data, 0.5/fs);
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
    preserve_percentage = 0.6;
    fc = 80;
    ITR = 150;
    T = size(data, 2);
    N_keep = round(T * preserve_percentage);
    kept_indexes = randi(T, [1, N_keep]);

    ecg_org = LPFilter(data(ch, :), fc / fs);
    ecg_dc = mean(ecg_org);
    ecg_zero_mean = ecg_org - ecg_dc;

    ecg_rec = zeros(1, T);
    ecg_rec(kept_indexes) = ecg_zero_mean(kept_indexes);

    ecg_all = zeros(ITR, T);
    for m = 1 : ITR
        ecg_all(m, :) = ecg_rec;
        ecg_rec(kept_indexes) = ecg_zero_mean(kept_indexes);
        ecg_rec_lp = LPFilter(ecg_rec, fc / fs);
        ecg_rec = ecg_rec_lp/sum(ecg_rec_lp.^2)*sum(ecg_rec.^2);
        %         hr_continuous = BPFilter5(hr_continuous, 10.0/fs, 35.0/fs, 3); % BPFilter5(x,fc,bw,order)
    end
    ecg_all = ecg_all + ecg_dc;
    ecg_rec = ecg_rec + ecg_dc;

    lgnd = {};
    t = (0 : T - 1) / fs;
    figure
    hold on
    plot(t, ecg_org, '-*', 'markersize', 14, 'linewidth', 2); lgnd = cat(2, lgnd, {'original'});
    %     plot(t, ecg_all');
    plot(t, ecg_all(end, :), 'k', 'linewidth', 2); lgnd = cat(2, lgnd, {'last iteration'});
    %     plot(t, ecg_all(1, :), 'b', 'linewidth', 2); lgnd = cat(2, lgnd, {'first iteration'});
    %     plot(t, (ecg_all(1, :) + ecg_all(end, :)) / 2.0, 'r', 'linewidth', 2);
    plot(t(kept_indexes), ecg_org(kept_indexes), 'ro', 'markersize', 16); lgnd = cat(2, lgnd, {'kept points'});
    grid
    legend(lgnd);


    %{
    f0 = 1.2;
    ch = 1;
    T = size(data, 2);
    peaks = PeakDetection(data(ch, :), f0/fs);
    peak_indexes = find(peaks);
    hr = 60 * fs ./ diff(peak_indexes);
    hr_indexes = peak_indexes(2:end);
    hr_dc = mean(hr);
    hr_zero_mean = hr - hr_dc;

%     hr_continuous = interp1([1, hr_indexes, T], [hr_zero_mean(1), hr_zero_mean, hr_zero_mean(end)], 1 : T, 'linear');

    hr_continuous = zeros(1, T); hr_continuous(hr_indexes) = hr_zero_mean;

    ITR = 50;
    fc = 45;
    HR = zeros(ITR, T);
    for m = 1 : ITR
        HR(m, :) = hr_continuous;
        hr_continuous(hr_indexes) = hr_zero_mean;
        hr_continuous_lp = LPFilter(hr_continuous, fc / fs);
        hr_continuous = hr_continuous_lp/sum(hr_continuous_lp.^2)*sum(hr_continuous.^2);
        %         hr_continuous = BPFilter5(hr_continuous, 10.0/fs, 35.0/fs, 3); % BPFilter5(x,fc,bw,order)
    end
    hr_continuous = hr_continuous + hr_dc;
    p = HR + hr_dc;

    % Plots results
    t = (0 : T - 1) / fs;
    figure
    plot(t, data(ch, :));
    hold on
    plot(t(peak_indexes), data(ch, peak_indexes), 'ro', 'markersize', 18);
    grid

    figure
    plot(t(hr_indexes), hr);
    grid

    figure
    %     plot(t, HR');
    hold on
    plot(t, HR(end, :), 'k', 'linewidth', 2);
    plot(t, HR(1, :), 'b', 'linewidth', 2);
    %     plot(t, (HR(1, :) + HR(end, :)) / 2.0, 'r', 'linewidth', 2);
    plot(t(hr_indexes), hr, 'ro', 'markersize', 16);
    grid
    %}
end


