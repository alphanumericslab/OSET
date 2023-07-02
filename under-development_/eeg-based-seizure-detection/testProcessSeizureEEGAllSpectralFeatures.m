close all;
clear;
clc;
resultfname = 'testProcessSeizureEEGAllSpectralFeatures3.txt';
% this is the same as the previous one; except that the feature vectors are
% kept in the frequency domain (look at F, F1, and F2 and how they are
% transposed). This gives equal number of features for all subjects
% regardless of the number of channels

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
ffs = 120;
f0 = 60; % powerline frequency Hz
Q = 60; % notch filter Q-factor
% wlens = [1.0 3.0 10.0]; % energy window lengths in seconds
% bins = 0:.01:3;
x = [];
drop = 1000; % the number of samples to exclude from the head and tail of the data during processing
% bins = 256;
nfft = 120;
SpectralWinLen = 5; % in seconds
SpectralOverlapPercentage = 50; % [0 99]
maxnumberofeigs = 15;

% preprocessing: resample the data to 120Hz with a polyphase bandpass filter passing [1Hz-55Hz]
load BandpassFilter1Hzto55HzAt400Hz h % load the bandpass filter used for resampling the input data
h1 = spectrum.welch('Hamming', round(SpectralWinLen*ffs), SpectralOverlapPercentage); % Create a Welch spectral estimator.

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
        
        Hx1 = MultiChannelSpectrum(x(:, 1:round(size(x,2)/2)), fs, nfft - 1, h1); % Estimate the PSD of the first half of the signal
        Hx2 = MultiChannelSpectrum(x(:, round(size(x,2)/2)+1:end), fs, nfft - 1, h1); % Estimate the PSD of the second part of the signal
        Hx = MultiChannelSpectrum(x, fs, nfft - 1, h1); % Estimate the PSD of the whole signal
        
        % % %         figure
        % % %         hold on
        % % %         for i = 1:size(x, 1),
        % % %             plot(Hx(i).Frequencies, log(Hx(i).data/sqrt(sum(Hx(i).data.^2))),'b');
        % % %         end
        % % %         grid;
        
        F1 = zeros(size(x,1), nfft/2);
        F2 = zeros(size(x,1), nfft/2);
        F = zeros(size(x,1), nfft/2);
        for i = 1:size(x, 1)
            F1(i, :) = (10*log10(Hx1(i).data/sqrt(sum(Hx1(i).data.^2)))).';
            F2(i, :) = (10*log10(Hx2(i).data/sqrt(sum(Hx2(i).data.^2)))).';
            F(i, :) = (10*log10(Hx(i).data/sqrt(sum(Hx(i).data.^2)))).';
        end
        
        %         Cf = F*F.'/size(F,1);
        %         Cfcross = F1*F2.'/size(F1,1);
        Cf = F.'*F/size(F,1);
        Cfcross = F1.'*F2/size(F1,1);
        
        feature = Cf(tril(true(size(Cf.')))); % take the lower triangular elements of the symmetric matrix
        feature_cross = Cfcross(:); % take all the elements of this matrix (it's not necessarily symmetric)
        egs = real(eig(Cf)); [~, II] = sort(egs, 1, 'descend');
        feature_eigs = 10*log10(egs(II)); % take the eigenvalues of the frequency domain correlation coefficients
        feature_eigs = feature_eigs(1:maxnumberofeigs); % take the first eigenvalues (up to the number of the channels)
        
        % write results
        fid = fopen(resultfname,'a');
        fprintf(fid, '%s\t%d\t%d\t%d', d(k).name, m, k, mode);
        for i = 1 : length(feature)
            fprintf(fid, '\t%8.4f', feature(i));
        end
        for i = 1 : length(feature_cross)
            fprintf(fid, '\t%8.4f', feature_cross(i));
        end
        for i = 1 : length(feature_eigs)
            fprintf(fid, '\t%8.4f', feature_eigs(i));
        end
        fprintf(fid, '\n');
        fclose(fid);
        disp(['subject: ', d(k).name]);
    end
end

clock