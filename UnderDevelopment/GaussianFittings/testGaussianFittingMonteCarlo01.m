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
            data = data - LPFilter(data, 1.5/fs);
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
    ITR = 100;
    samples = 1 : round(1.0*fs);
    sd = 0.01 * std(data(ch, :));
    A = zeros(ITR, size(data(ch, :), 2));
    for itr = 1 : ITR
        x = data(ch, :) + sd * randn(size(data(ch, :)));
        A(itr, :) = diff([x(1) x]) ./ x;
    end
    A_mn = mean(A, 1);
    A_std = std(A, [], 1);
    figure
    subplot(211)
    plot(data(ch, samples));
%     plot(A', 'r');
    grid
    subplot(212)
    errorbar(A_mn(samples), A_std(samples)/2);
    hold on
    plot(A_mn(samples));
    grid
    
    color_code = tanh(A_std/(0.5*median(A_std)));
    figure
    scatter(x, diff([x(1) x]), 10, (ones(3, 1)*color_code)', 'filled');
    grid
end

