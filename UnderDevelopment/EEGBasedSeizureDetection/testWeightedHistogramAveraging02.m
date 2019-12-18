close all;
clear;
clc;
resultfname = 'testWeightedHistogramAveraging02.txt';

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
f0 = 60; % powerline frequency Hz
Q = 30; % notch filter Q-factor
wlen = [0.5 1.0 3.0 5.0 10.0]; % energy window lengths in seconds
drop = 200; % remove data head and tail due to filtering transitions
ffs = 160; % resample rate

% bins = 0:.01:3;
% epochlen = 15; % epoch length in seconds
% epochoverlap = 12; % overlap between epochs in seconds
% smoothing_iterations = 100;
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
        %         PlotECG(x, 4, 'b', fs);
        
        f = zeros(1, length(wlen));
        for i = 1:length(wlen)
            ll = round(wlen(i)*fs);
            E = sqrt(filter(ones(1, ll), ll, x.^2, [], 2));
            
            % remove data head and tail due to filtering transitions
            E = E(:, drop:end-drop);
            
            % normalize
            EE = (E - mean(E, 2)*ones(1, size(E,2)))./(std(E, [], 2)*ones(1, size(E,2)));
            
            Cx = cov(EE');
            I = find(tril(Cx, -1));
            f(i) = mean(Cx(I));
        end
        
        % write results
        fid = fopen(resultfname,'a');
        fprintf(fid, '%s\t%d\t%d\t%d\t', d(k).name, m, k, mode);
        for i = 1:length(wlen)
            fprintf(fid, '%8.6f\t', f(i));
        end
        fprintf(fid, '\n');
        fclose(fid);
        disp(['subject: ', d(k).name]);
    end
end

% for i = 1 : size(x, 1)
%     figure;
%     t = (0:length(x)-1)/fs;
%     plot(t, x(i, :),'k');
%     hold on
%     t = (0:length(E)-1)/fs;
%     plot(t, E(i, :));
%     grid
% end
% 
% figure
% imshow(Cx);
% grid
% 
clock
