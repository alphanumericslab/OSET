close all;
clear;
clc;
resultfname = 'testProcessSeizureEEGICAOnNormalizedHistogram.txt';

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
ffs = 120;
f0 = 60; % powerline frequency Hz
Q = 60; % notch filter Q-factor
% wlens = [1.0 3.0 10.0]; % energy window lengths in seconds
% bins = 0:.01:3;
x = [];
drop = 500; % the number of samples to exclude from the head and tail of the data during processing
bins = 256;

% preprocessing: resample the data to 120Hz with a polyphase bandpass filter passing [1Hz-55Hz]
load BandpassFilter1Hzto55HzAt400Hz h % load the bandpass filter used for resampling the input data

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
        if(fs < 500)
            x = resample(x', ffs, round(fs), h*round(fs)/ffs)';
        else
            midfs = 400;
            % resample to intermediate frequency
            x = resample(x', midfs, round(fs))';
            
            % IIR notch
            Wo = f0/(midfs/2);  BW = Wo/Q;
            [b,a] = iirnotch(Wo, BW);
            x = filter(b, a, x, [], 2);
            
            % resample using bandpass filter
            x = resample(x', ffs, midfs, h*midfs/ffs)';
        end
        fs = ffs;
        
        % drop transient samples
        x = x(:,drop:end-drop);
        
        % normalize
        x = (x - mean(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
        % x = (x - median(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
        
        Hx = zeros(size(x,1), bins);
        Hvalsx = zeros(size(x,1), bins);
        for i = 1:size(x, 1)
            [Hx(i, :) Hvalsx(i, :)] = hist(x(i, :), bins);
        end
        
        cdfx = cumsum(Hx, 2); cdfx = cdfx./(cdfx(:,end)*ones(1, size(cdfx,2)));
        
        % xx is x with a normalized histogram
        xx = zeros(size(x));
        for i = 1:size(x, 1)
            for j = 1:size(x, 2)
                I = find(Hvalsx(i, :) >= x(i, j), 1, 'first');
                if(isempty(I))
                    I = size(cdfx,2);
                end
                xx(i, j) = cdfx(i, I) - 0.5;
            end
        end
        wxJADE = jadeR(xx);
        sx = wxJADE*xx;
        
        Hsx1 = zeros(size(sx,1), bins);
        Hsx2 = zeros(size(sx,1), bins);
        Hsvalsx1 = zeros(size(sx,1), bins);
        Hsvalsx2 = zeros(size(sx,1), bins);
        erx = zeros(1, size(sx,1));
        erxn = zeros(1, size(sx,1));
        mn = zeros(1, size(sx,1));
        md = zeros(1, size(sx,1));
        sd = zeros(1, size(sx,1));
        skw = zeros(1, size(sx,1));
        krt = zeros(1, size(sx,1));
        % find assymetry in the histogram of the first and second half of sx and sy
        for i = 1:size(sx, 1)
            [Hsx1(i, :) Hsvalsx1(i, :)] = hist(sx(i, 1:round(size(sx,2)/2)), bins);
            Hsx2(i, :) = hist(sx(i, round(size(sx,2)/2)+1:end), Hsvalsx1(i, :));
            erx(i) = mean(abs(Hsx1(i, :) - Hsx2(i, :)));
            erxn(i) = mean(abs(Hsx1(i, :)/sum(Hsx1(i, :)) - Hsx2(i, :)/sum(Hsx2(i, :))));
            mn(i) = mean(sx(i,:));
            md(i) = median(sx(i,:));
            sd(i) = std(sx(i,:), 1);
            skw(i) = skewness(sx(i,:));
            krt(i) = kurtosis(sx(i,:)) - 3.0;
        end
        errorx = sort(erx(:), 1, 'descend');
        errorxn = sort(erxn(:), 1, 'descend');
        mnx = sort(mn(:), 1, 'descend');
        mdx = sort(md(:), 1, 'descend');
        sdx = sort(sd(:), 1, 'descend');
        skwx = sort(skw(:), 1, 'descend');
        krtx = sort(krt(:), 1, 'descend');
        
        f1 = mean(errorx(1:4));
        f2 = mean(errorxn(1:4));
        f3 = mean(mnx(1:4));
        f4 = mean(mdx(1:4));
        f5 = mean(sdx(1:4));
        f6 = mean(skwx(1:4));
        f7 = mean(krtx(1:4));
        
        f8 = median(errorx);
        f9 = median(errorxn);
        f10 = median(mnx);
        f11 = median(mdx);
        f12 = median(sdx);
        f13 = median(skwx);
        f14 = median(krtx);

        % write results
        fid = fopen(resultfname,'a');
        fprintf(fid, '%s\t\t%d\t%d\t%d\t', d(k).name, m, k, mode);
        fprintf(fid, '%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t', f1, 1e4*f2, f3, f4, f5, f6, f7);
        fprintf(fid, '%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n', f8, 1e4*f9, f10, f11, f12, f13, f14);
        fclose(fid);
        disp(['subject: ', d(k).name]);
    end
end

clock