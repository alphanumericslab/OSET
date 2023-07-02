close all;
clear;
clc;
resultfname = 'testProcessSeizureEEG_EnergyDistributions02.txt';

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
ffs = 150; % resample rate
f0 = 60; % powerline frequency Hz
Q = 30; % notch filter Q-factor
wlens = [1.0 3.0 10.0]; % energy window lengths in seconds
bins = 0:.01:3;
x = [];
for m = 1 : length(subject),
    % search path and subjects
    path = ['E:\Sameni\Projects\Seizure\' subject{m} '\'];
    d = dir([path '*.mat']);
    NumRecords = length(d);
    for k = 1 : NumRecords,
        % clear previous variable to save memory
        if(~isempty(x))
            clear x;
        end
        
        % load data
        load([path d(k).name]);
        reg = regexp(d(k).name, '_');
        varname = [d(k).name(reg(2)+1:reg(4)-1) '_' num2str(str2double(d(k).name(end-7:end-4)))];
        md = d(k).name(reg(2)+1:reg(3)-1);
        if(isequal(md, 'interictal'))
            mode = 1;
        elseif(isequal(md, 'preictal'))
            mode = 2;
        elseif(isequal(md, 'test'))
            mode = 3;
        end
        eval(['data = ' varname '; clear ' varname ';']);
        x = data.data;
        fs = data.sampling_frequency;
        clear data;
        
        % resample data
        if(fs > ffs) % for the patients who have high sampling rates
            x = resample(x', ffs, round(fs))';
            fs = ffs;
        end
        
        % pre-process
        % x = LPFilter(x, 80.0/fs); % lowpass filter
        x = x - LPFilter(x, 0.1/fs); % highpass filter
        % IIR notch
        Wo = f0/(fs/2);  BW = Wo/Q;
        [b,a] = iirnotch(Wo,BW);
        x = filter(b, a, x, [], 2);
        
        % normalize
        x = (x - mean(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
        
        % remove data head and tail due to transitions
        drop = 200;
        xx = x(:, drop:end-drop);
        
        channelnum = size(xx,1);
        WinN = length(wlens);
        E = zeros(WinN, channelnum, size(xx,2));
        histogram = zeros(WinN, length(bins));
        
        mn = zeros(1, WinN);
        md = zeros(1, WinN);
        sd = zeros(1, WinN);
        skw = zeros(1, WinN);
        krt = zeros(1, WinN);
        for i = 1:WinN
            ll = round(wlens(i)*fs);
            E(i, :, :) = sqrt(filter(ones(1, ll), ll, xx.^2, [], 2));
            e = squeeze(E(i, :, :));
            [nn, ee] = hist(e(:), bins);
            NN = sum(nn);
            histogram(i, :) = nn/NN;
            
            mn(i) = mean(histogram(i, :));
            md(i) = median(histogram(i, :));
            sd(i) = std(histogram(i, :));
            skw(i) = skewness(histogram(i, :));
            krt(i) = kurtosis(histogram(i, :)) - 3;
            
%             figure;
%             plot(ee, nn/NN, 'color', (i-1)*ones(1,3)/WinN);
%             grid;
        end
        % write results
        fid = fopen(resultfname,'a');
        fprintf(fid, '%s\t%d\t%d\t%d\t', d(k).name, m, k, mode);
        for i = 1:WinN
            fprintf(fid, '%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t', mn(i), md(i), sd(i), skw(i), krt(i));
        end
        fprintf(fid, '\n');
        fclose(fid);
        disp(['subject: ', d(k).name]);
    end
end

for i = 1 : size(x, 1)
    figure;
    t = (0:length(x)-1)/fs;
    plot(t, x(i, :),'k');
    hold on
    t = (0:length(E)-1)/fs;
    plot(t, squeeze(E(:, i, :)));
    grid
end

clock