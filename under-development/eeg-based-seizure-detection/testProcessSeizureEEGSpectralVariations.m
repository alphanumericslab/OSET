close all;
clear;
clc;

rawspectra = 'AllSubjectsRawSpectra.txt';
resultfname1 = 'testProcessSeizureEEGRWAofRawSpectra.txt';
resultfname2 = 'testProcessSeizureEEGAllSpectralFeatures5.txt';

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
ffs = 120;
f0 = 60; % powerline frequency Hz
Q = 60; % notch filter Q-factor
x = [];
drop = 1000; % the number of samples to exclude from the head and tail of the data during processing
nfft = 120;
SpectralWinLen = 10; % in seconds
SpectralOverlapPercentage = 75; % [0 99]
per_subject = 1;

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
        
        %         Hx1 = MultiChannelSpectrum(x(:, 1:round(size(x,2)/2)), fs, nfft - 1, h1); % Estimate the PSD of the first half of the signal
        %         Hx2 = MultiChannelSpectrum(x(:, round(size(x,2)/2)+1:end), fs, nfft - 1, h1); % Estimate the PSD of the second part of the signal
        Hx = MultiChannelSpectrum(x, fs, nfft - 1, h1); % Estimate the PSD of the whole signal
        
        %         close all
        % write results
        fidraw = fopen(rawspectra,'a');
%         figure
%         hold on
        HHH = zeros(size(x, 1), nfft/2);
        for i = 1:size(x, 1),
            fprintf(fidraw, '%s\t%d\t%d\t%d\t%d', d(k).name, m, k, mode, i);
            %             fprintf(fidraw, '\t%10.6f', Hx1(i).data/sum(Hx1(i).data.^2));
            %             fprintf(fidraw, '\t%10.6f', Hx2(i).data/sum(Hx2(i).data.^2));
            HH = Hx(i).data/sum(Hx(i).data.^2);
            for j = 1:length(HH),
                fprintf(fidraw, '\t%10.6f', HH(j));
            end
            fprintf(fidraw, '\n');
            
            HHH(i, :) = HH';
%             plot(Hx(i).Frequencies, log(HH),'b');
        end
        fclose(fidraw);
        HHHAVG = RWAverage(HHH);
%         plot(Hx(1).Frequencies, log(HHHAVG),'r', 'linewidth', 3);
%         grid
        
        fidrwa = fopen(resultfname1,'a');
        fprintf(fidrwa, '%s\t%d\t%d\t%d', d(k).name, m, k, mode);
        for i = 1:length(HHHAVG),
            fprintf(fidrwa, '\t%10.6f', log(HHHAVG(i)));
        end
        fprintf(fidrwa, '\n');
        fclose(fidrwa);

        disp(['subject: ', d(k).name]);
    end
end

% % % spc = importdata(rawspectra);
% % %
% % % filenames_all = spc.textdata;
% % % subjects_all = spc.data(:, 1);
% % % trials_all = spc.data(:, 2);
% % % types_all = spc.data(:, 3);
% % % chanels_all = spc.data(:, 4);
% % % % spectra1 = spc.data(:, 5:4+nfft/2);
% % % % spectra2 = spc.data(:, 5+nfft/2:4+nfft);
% % % % spectra = spc.data(:, 5+nfft:4+3*nfft/2);
% % % spectra = spc.data(:, 5:end);
% % %
% % % interictal_label = 1;
% % % preictal_label = 2;
% % % test_label = 3;
% % %
% % % if(per_subject == 1)
% % %     max_subject = max(subjects_all);
% % % else
% % %     max_subject = 1;
% % % end
% % %
% % % % scores = [];
% % % for s = 1 : max_subject,
% % %     if(per_subject == 1)
% % %         subjects = (subjects_all == s);
% % %     else
% % %         subjects = true(size(subjects_all));
% % %     end
% % %
% % %     types = types_all(subjects); % types (interictal/preictal/test) per subject
% % %     trials = trials_all(subjects);
% % %     channels = chanels_all(subjects);
% % %
% % %     %     subjectSpectra1 = spectra1(subjects, :);
% % %     %     subjectSpectra2 = spectra2(subjects, :);
% % %     subjectSpectra = spectra(subjects, :);
% % %
% % %     interictal_indexes = find(types == interictal_label); % interictal
% % %     preictal_indexes = find(types == preictal_label); % preictal
% % %     train_indexes = [interictal_indexes ; preictal_indexes];
% % %     test_indexes = find(types == test_label); % test
% % %
% % %     %     AverageSpectra1Interictal = zeros(max(channels), nfft/2);
% % %     %     AverageSpectra2Interictal = zeros(max(channels), nfft/2);
% % %     AverageSpectraInterictal = zeros(max(channels), nfft/2);
% % %
% % %     %     AverageSpectra1Preictal = zeros(max(channels), nfft/2);
% % %     %     AverageSpectra2Preictal = zeros(max(channels), nfft/2);
% % %     AverageSpectraPreictal = zeros(max(channels), nfft/2);
% % %     for i = 1 : max(channels),
% % %         %         AverageSpectra1Interictal(i, :) = RWAverage(subjectSpectra1(channels == i && types == interictal_label, :), :);
% % %         %         AverageSpectra2Interictal(i, :) = RWAverage(subjectSpectra2(channels == i && types == interictal_label, :), :);
% % %         AverageSpectraInterictal(i, :) = RWAverage(subjectSpectra(channels == i && types == interictal_label, :), :);
% % %
% % %         %         AverageSpectra1Preictal(i, :) = RWAverage(subjectSpectra1(channels == i && types == preictal_label, :), :);
% % %         %         AverageSpectra2Preictal(i, :) = RWAverage(subjectSpectra2(channels == i && types == preictal_label, :), :);
% % %         AverageSpectraPreictal(i, :) = RWAverage(subjectSpectra(channels == i && types == preictal_label, :), :);
% % %     end
% % %
% % %     NumRecords = max(trials);
% % %     for k = 1 : NumRecords,
% % %         for i = 1 : max(channels),
% % %             df = subjectSpectra(channels == i && trials == k, :) - AverageSpectraInterictal(i, :);
% % %             er_mean(k) = mean()
% % %         end
% % %     end
% % %
% % %
% % %     % % %         figure
% % %     % % %         hold on
% % %     % % %         for i = 1:size(x, 1),
% % %     % % %             plot(Hx(i).Frequencies, log(Hx(i).data/sqrt(sum(Hx(i).data.^2))),'b');
% % %     % % %         end
% % %     % % %         grid;
% % %
% % %     F1 = zeros(size(x,1), nfft/2);
% % %     F2 = zeros(size(x,1), nfft/2);
% % %     F = zeros(size(x,1), nfft/2);
% % %     for i = 1:size(x, 1)
% % %         F1(i, :) = (10*log10(Hx1(i).data/sqrt(sum(Hx1(i).data.^2)))).';
% % %         F2(i, :) = (10*log10(Hx2(i).data/sqrt(sum(Hx2(i).data.^2)))).';
% % %         F(i, :) = (10*log10(Hx(i).data/sqrt(sum(Hx(i).data.^2)))).';
% % %     end
% % %
% % %     Cf = F.'*F/size(F,1);
% % %     Cfcross = F1.'*F2/size(F1,1);
% % %
% % %     feature = Cf(tril(true(size(Cf.')))); % take the lower triangular elements of the symmetric matrix
% % %     feature_cross = Cfcross(:); % take all the elements of this matrix (it's not necessarily symmetric)
% % %     egs = real(eig(Cf)); [~, II] = sort(egs, 1, 'descend');
% % %     feature_eigs = 10*log10(egs(II)); % take the eigenvalues of the frequency domain correlation coefficients
% % %     feature_eigs = feature_eigs(1:maxnumberofeigs); % take the first eigenvalues (up to the number of the channels)
% % %
% % %     L = 4;
% % %     wtype = {'db4','sym5','coif3'};
% % %     waveletfeaturesCf = [];
% % %     for i = 1:length(wtype),
% % %         [C, S] = wavedec2(Cf, L, wtype{i});
% % %         for LL = 1:L;
% % %             A = appcoef2(C, S, wtype{i}, LL);
% % %             waveletfeaturesCf = [waveletfeaturesCf log(sqrt(mean(A(:).^2)))];
% % %         end
% % %         [Ea, Eh, Ev, Ed] = wenergy2(C, S);
% % %         waveletfeaturesCf = [waveletfeaturesCf log([Ea(:) ; Eh(:) ; Ev(:) ; Ed(:)]')];
% % %     end
% % %
% % %     waveletfeaturesCfcross = [];
% % %     for i = 1:length(wtype),
% % %         [C, S] = wavedec2(Cfcross, L, wtype{i});
% % %         for LL = 1:L;
% % %             A = appcoef2(C, S, wtype{i}, LL);
% % %             waveletfeaturesCfcross = [waveletfeaturesCfcross log(sqrt(mean(A(:).^2)))];
% % %         end
% % %         [Ea, Eh, Ev, Ed] = wenergy2(C, S);
% % %         waveletfeaturesCfcross = [waveletfeaturesCfcross log([Ea(:) ; Eh(:) ; Ev(:) ; Ed(:)]')];
% % %     end
% % %
% % %     % write results
% % %     fid = fopen(resultfname,'a');
% % %     fprintf(fid, '%s\t%d\t%d\t%d', d(k).name, m, k, mode);
% % %     for i = 1 : length(feature)
% % %         fprintf(fid, '\t%10.6f', feature(i));
% % %     end
% % %     for i = 1 : length(feature_cross)
% % %         fprintf(fid, '\t%10.6f', feature_cross(i));
% % %     end
% % %     for i = 1 : length(feature_eigs)
% % %         fprintf(fid, '\t%10.6f', feature_eigs(i));
% % %     end
% % %     fprintf(fid, '\n');
% % %     fclose(fid);
% % %
% % %     fid_wavelet = fopen(resultfname_wavelet,'a');
% % %     fprintf(fid_wavelet, '%s\t%d\t%d\t%d', d(k).name, m, k, mode);
% % %     for i = 1 : length(waveletfeaturesCf)
% % %         fprintf(fid_wavelet, '\t%10.6f', waveletfeaturesCf(i));
% % %     end
% % %     for i = 1 : length(waveletfeaturesCfcross)
% % %         fprintf(fid_wavelet, '\t%10.6f', waveletfeaturesCfcross(i));
% % %     end
% % %     fprintf(fid_wavelet, '\n');
% % %     fclose(fid_wavelet);
% % %
% % %     disp(['subject: ', d(k).name]);
% % % end
% % % end

clock