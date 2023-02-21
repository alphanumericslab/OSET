% Test script for SpatioTemporalInnovationsFilterDesigner, which estimates an innovations filter from ensembles of multichannel random
% processes.
%
% The Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET
% Reza Sameni, Feb 2023

close all
clear
clc

case_study = 'EEG'; % 'EEG' or 'PCG'
switch case_study
    case 'PCG' % The CirCor PCG dataset on PhysioNet
        datafilepath = '../../../DataFiles/physionet.org/files/circor-heart-sound/1.0.3/'; % Change this path to where you have the .wav data files
        filelist = dir(fullfile([datafilepath, '**/*.wav']));  % get list of all mat files of interest
        %     filelist = filelist(1:100);

        sample_data = cell(1, length(filelist));
        for k = 1 : length(filelist) % Sweep over all or specific records
            datafilename = [filelist(k).folder '/' filelist(k).name];
            headerfilename = [datafilename(1:end-3) 'hea'];
            fid = fopen(headerfilename);
            tline = fgetl(fid);
            fclose(fid);
            tline_split = split(tline, ' ');
            fs = str2double(tline_split{3});

            [record, Fs] = audioread(datafilename);
            sample_data{k} = record';
            disp(k)
        end

    case 'EEG' % A sample EEG from OSET
        load EEGdata.mat
        fs = 250; % Sampling frequency
        num_segments = 5; % Number of segments used to split the large data matrix
        wlen = round(fs * 10.0); % length of each segment
        sample_data = cell(1, num_segments);
        for k = 1 : num_segments
            sample_data{k} = data((k-1)*wlen + 1 : k*wlen, :)';
        end
    otherwise
        error('Undefined case study. Try EEG, PCG, or adapt the script according to your data.')
end

params.spatial_filter_type = 'BY_PASS'; % 'BY_PASS', 'PCA' or 'ICA'
params.normalize_records = true; % true/false
params.fs = fs; % sampling frequency
params.keep_mean = true; % true/false
params.spectral_len = 513; % number of frequency bins for spectral estimation
params.filter_len = params.spectral_len; % innovations filter length (best practice to set equal to params.spectral_len)
params.smooth_spectrum = false; % true/false smooth the spectrum before spectral factorization
if params.smooth_spectrum
    params.lambda = 10000.0; % Tikhonov regularization factor used for spectral smoothing
end
params.spectral_averaging_method = 'MEDIAN'; % 'MEDIAN', 'MEAN', 'MAX', 'MIN', 'MAX_MIN_AVG'
params.innovation_filter_type = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE';
params.plot_results = true; % true/false plot results

% Design the innovations filter
[h_innovations, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = SpatioTemporalInnovationsFilterDesigner(sample_data, params);


% Test the innovations filter
dat = sample_data{1};
syn_signal_len = size(dat, 2);
x = randn(length(h_innovations), syn_signal_len);
y = randn(length(h_innovations), syn_signal_len);
for ch = 1 : length(h_innovations)
    y(ch, :) = filter(h_innovations{ch}, 1, x(ch, :));
    [H,F] = freqz(h_innovations{ch}, 1, [], fs);

    figure
    subplot(411)
    plot((0 : length(dat(ch, :)) - 1)/fs, dat(ch, :));
    grid
    xlabel('time(s)');
    ylabel('real data');
    set(gca, 'fontsize', 18)

    subplot(412)
    plot((0 : length(y(ch, :)) - 1)/fs, y(ch, :));
    grid
    xlabel('time(s)');
    ylabel('innovations process');
    set(gca, 'fontsize', 18)

    subplot(413)
    plot(F, 20*log10(abs(H)));
    grid
    xlabel('frequency(Hz)');
    ylabel('innovations filter magnitude');
    set(gca, 'fontsize', 18)

    subplot(414)
    plot(F, unwrap(angle(H)));
    grid
    xlabel('frequency(Hz)');
    ylabel('innovations filter phase');
    set(gca, 'fontsize', 18)
    sgtitle('Real data vs synthetic innovations process', 'fontsize', 18);
end

% save InnovationsFilterDesign h
