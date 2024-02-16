% A test script to extract R-peaks from sample records of the PTB-XL
% dataset

close all
clear
clc

data_path = '../sample-data/ptb-xl/*.mat'; % set data path
fs = 500; % set sampling frequency

data_list = dir(data_path);
for k = 1 : length(data_list)
    fname = fullfile(data_list(k).folder, data_list(k).name);
    data = load(fname).data;

    % R-peak detector
    % parameters may need tweaking for improved performance per channel
    params = [];
    params.min_peak_distance = 0.3;
    params.p_residual_th = 95.0;
    params.p_signal_th = 90.0;
    [~, I_peaks] = peak_det_likelihood(data, fs, params);

    % heart rate
    hr = 60*fs./diff(I_peaks);

    % plot results
    t = (0:length(data)-1)/fs;

    figure
    hold on
    plot(t, data)
    plot(t(I_peaks), data(I_peaks), 'o', 'MarkerFaceColor', 'g', 'markersize', 16)
    grid
    xlabel('time(s)')
    ylabel('Amplitude(normalized)')
    set(gca, 'fontsize', 18)
    title('ECG and R-peaks');

    figure
    plot(t(I_peaks(2:end)), hr)
    grid
    xlabel('time(s)')
    ylabel('Heart rate (bpm)')
    set(gca, 'fontsize', 18)
    title('Heart rate');

end

