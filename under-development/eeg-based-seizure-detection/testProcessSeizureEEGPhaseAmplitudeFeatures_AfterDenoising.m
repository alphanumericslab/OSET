close all;
clear;
clc;
resultfname_amp = 'testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_amp.txt';
resultfname_phase = 'testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_phase.txt';
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
% h1 = spectrum.welch('Hamming', round(SpectralWinLen*ffs), SpectralOverlapPercentage); % Create a Welch spectral estimator.

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
        
        % remove spikes and nonstationarities from the data
        % % %         x = RemoveEEGSeizureNoise(x, fs);
        y = RemoveEEGSeizureNoise(x, fs, 0.8, 2.0, 10);
        
        X = fft(x, nfft, 2);
        X1 = fft(x(:, 1:round(size(x,2)/2)), nfft, 2);
        X2 = fft(x(:, round(size(x,2)/2)+1:end), nfft, 2);
        
        Xamp = abs(X);
        X1amp = abs(X1);
        X2amp = abs(X2);
        
        Xphase = zeros(size(X,1),nfft/2);
        X1phase = zeros(size(X1,1),nfft/2);
        X2phase = zeros(size(X2,1),nfft/2);
        for i = 1:size(x, 1),
            Xphase(i, :) = phase(X(i, 1:nfft/2)) - LPFilter(phase(X(i, 1:nfft/2)), 1.0/fs);
            X1phase(i, :) = phase(X1(i, 1:nfft/2)) - LPFilter(phase(X1(i, 1:nfft/2)), 1.0/fs);
            X2phase(i, :) = phase(X2(i, 1:nfft/2)) - LPFilter(phase(X2(i, 1:nfft/2)), 1.0/fs);
        end
        freq = fs*(0:nfft-1)/nfft;
        
        Famp = 20*log10(Xamp(:, 1:nfft/2));
        F1amp = 20*log10(X1amp(:, 1:nfft/2));
        F2amp = 20*log10(X2amp(:, 1:nfft/2));
        
        C_amp = Famp.'*Famp/size(Famp,1);
        CCross_amp = F1amp.'*F2amp/size(F1amp,1);
        
        C_phase = Xphase.'*Xphase/size(Xphase,1);
        CCross_phase = X1phase.'*X2phase/size(X1phase,1);
        
% % %                 figure
% % %                 subplot(121);
% % %                 mesh(C_amp);
% % % % % %                 hold on
% % % % % %                 for i = 1:size(x, 1),
% % % % % %                     plot(freq, 20*log10(Xamp),'b');
% % % % % %                     plot(freq, 20*log10(X1amp),'r');
% % % % % %                     plot(freq, 20*log10(X2amp),'g');
% % % % % %                 end
% % % % % %                 grid;
% % %                 subplot(122);
% % %                 mesh(C_phase);
% % % % % %                 hold on
% % % % % %                 for i = 1:size(x, 1),
% % % % % %                     plot(freq, Xphase,'b');
% % % % % %                     plot(freq, X1phase,'r');
% % % % % %                     plot(freq, X2phase,'g');
% % % % % %                 end
% % % % % %                 grid;
        
        feature_amp = C_amp(tril(true(size(C_amp.')))); % take the lower triangular elements of the symmetric matrix
        feature_cross_amp = CCross_amp(:); % take all the elements of this matrix (it's not necessarily symmetric)
        egs_amp = real(eig(C_amp)); [~, II] = sort(egs_amp, 1, 'descend');
        feature_eigs_amp = 10*log10(egs_amp(II)); % take the eigenvalues of the frequency domain correlation coefficients
        feature_eigs_amp = feature_eigs_amp(1:maxnumberofeigs); % take the first eigenvalues (up to the number of the channels)
        
        feature_phase = C_phase(tril(true(size(C_phase.')))); % take the lower triangular elements of the symmetric matrix
        feature_cross_phase = CCross_phase(:); % take all the elements of this matrix (it's not necessarily symmetric)
        egs_phase = real(eig(C_phase)); [~, II] = sort(egs_phase, 1, 'descend');
        feature_eigs_phase = 10*log10(egs_phase(II)); % take the eigenvalues of the frequency domain correlation coefficients
        feature_eigs_phase = feature_eigs_phase(1:maxnumberofeigs); % take the first eigenvalues (up to the number of the channels)
        
        % write results
        fid_amp = fopen(resultfname_amp,'a');
        fid_phase = fopen(resultfname_phase,'a');
        fprintf(fid_amp, '%s\t%d\t%d\t%d', d(k).name, m, k, mode);
        fprintf(fid_phase, '%s\t%d\t%d\t%d', d(k).name, m, k, mode);
        for i = 1 : length(feature_amp)
            fprintf(fid_amp, '\t%8.4f', feature_amp(i));
            fprintf(fid_phase, '\t%8.4f', feature_phase(i));
        end
        for i = 1 : length(feature_cross_amp)
            fprintf(fid_amp, '\t%8.4f', feature_cross_amp(i));
            fprintf(fid_phase, '\t%8.4f', feature_cross_phase(i));
        end
        for i = 1 : length(feature_eigs_amp)
            fprintf(fid_amp, '\t%8.4f', feature_eigs_amp(i));
            fprintf(fid_phase, '\t%8.4f', feature_eigs_phase(i));
        end
        fprintf(fid_amp, '\n');
        fprintf(fid_phase, '\n');
        fclose(fid_amp);
        fclose(fid_phase);
        disp(['subject: ', d(k).name]);
    end
end

clock