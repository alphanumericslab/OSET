close all;
clear;
clc;
resultfname = 'results.txt';

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
f0 = 60; % powerline frequency Hz
Q = 30; % notch filter Q-factor
wlen = 5.0; % energy window lengths in seconds
bins = 0:.01:3;
epochlen = 15; % epoch length in seconds
epochoverlap = 12; % overlap between epochs in seconds
smoothing_iterations = 100;
x = [];
for m = 4 : 4,%length(subject),
    % search path and subjects
    path = ['E:\Sameni\Projects\Seizure\' subject{m} '\'];
    d = dir([path '*.mat']);
    NumRecords = length(d);
    for k = 700 : 700,%NumRecords,
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
        ffs = 250; % resample rate
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
        x = x(:, drop:end-drop);
        %         PlotECG(x, 4, 'b', fs);
        
        ll = round(wlen*fs);
        %         E = abs(x);
        E = sqrt(filter(ones(1, ll), ll, x.^2, [], 2));
        
        % normalize
        EE = (E - mean(E, 2)*ones(1, size(E,2)))./(std(E, [], 2)*ones(1, size(E,2)));
        Cx = cov(EE');
        % Cx = E*E'/size(E, 2);
        
        channelnum = size(x,1);
        epochnum = floor((length(x)-epochlen*fs)/((epochlen-epochoverlap)*fs));
        
        for p = 1:channelnum,
            histogram = zeros(length(bins), epochnum);
            for i = 1:epochnum,
                e = E(:, round(i*(epochlen-epochoverlap)*fs:i*(epochlen-epochoverlap)*fs+epochlen*fs));
                [nn, ee] = hist(e(p, :), bins);
                NN = sum(nn);
                histogram(:, i) = nn/NN;
            end
            
            w = ones(1, epochnum)/epochnum;
            figure;
            hold on
            for s = 1 : smoothing_iterations, % smoothing iterations
                meanhist = sum(histogram.*w(ones(1,length(bins)), :), 2);
                error = histogram - meanhist(:, ones(1, size(histogram, 2)));
                dev = std(error, [], 1);
                %                 w = dev - min(dev);
                %                 w = 1 - w./sum(w);
                %                 w = dev./sum(dev);
                w = 1./dev;
                w = w ./ sum(w);
                plot(bins, meanhist);%, 'k', 'linewidth', 2);
            end
            title(num2str(p));
            grid
        end
        %         [Y,I] = sort(dev, 1, 'ascend');
        %             figure;
        %             plot(ee, nn/NN, 'color', (i-1)*ones(1,3)/WinN);
        %             grid;
        % write results
        %         fid = fopen(resultfname,'a');
        %         fprintf(fid, '%s\t%d\t%d\t%d\t', d(k).name, m, k, mode);
        %         for i = 1:WinN
        %             fprintf(fid, '%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t', mn(i), md(i), sd(i), skw(i), krt(i));
        %         end
        %         fprintf(fid, '\n');
        %         fclose(fid);
        disp(['subject: ', d(k).name]);
    end
end

figure
imshow(Cx);
grid

for i = 1 : size(x, 1)
    figure;
    t = (0:length(x)-1)/fs;
    plot(t, x(i, :),'k');
    hold on
    t = (0:length(E)-1)/fs;
    plot(t, E(i, :));
    grid
end

% clock