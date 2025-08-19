clear
close all
clc

spectral_window = 5.0;
spectral_overlap = 4.75;
spectral_len = 1000;
spectral_overlap_len = round(spectral_len*0.5);

dataset1_path = '/Users/rsameni/Downloads/013000/';
dataset2_path = '/Users/rsameni/Downloads/g1/';

dataset1_fnames = dir(fullfile(dataset1_path, '*.mat'));
dataset2_fnames = dir(fullfile(dataset2_path, '*.mat'));

% ch = 3;

% Sx = zeros(length(dataset1_fnames), spectral_len);
Sx = [];
for k = 1:100%length(dataset1_fnames)
    record = fullfile(dataset1_fnames(k).folder, dataset1_fnames(k).name(1:end-4));
    record_desc = wfdb_desc(record);
    for ch = 1 : length(record_desc)
        fs = record_desc(ch).SamplingFrequency;
        data = load(record).val;
        if length(data(ch, :)) >= spectral_len && fs == 500.0
            % data(ch, :) = data(ch, :) - mean(data(ch, :));
            data(ch, :) = data(ch, :) - lp_filter_zero_phase(data(ch, :), 0.5/fs);
            % Sx(k, :) = pwelch(data(ch, :), hamming(spectral_len), spectral_overlap_len, spectral_len, fs, 'twosided')';
            Sx = cat(1, Sx, pwelch(data(ch, :), hamming(spectral_len), spectral_overlap_len, spectral_len, fs, 'twosided')');
        end
    end
end
Sx = median(Sx, 1);
Sx = Sx / sum(Sx);

% Sy = zeros(length(dataset2_fnames), spectral_len);
Sy = [];
for k = 1:100%length(dataset2_fnames)
    record = fullfile(dataset2_fnames(k).folder, dataset2_fnames(k).name(1:end-4));
    record_desc = wfdb_desc(record);
    for ch = 1 : length(record_desc)
        fs = record_desc(ch).SamplingFrequency;
        data = load(record).val;
        if length(data(ch, :)) >= spectral_len  && fs == 500.0
            % data(ch, :) = data(ch, :) - mean(data(ch, :));
            data(ch, :) = data(ch, :) - lp_filter_zero_phase(data(ch, :), 0.5/fs);
            % Sy(k, :) = pwelch(data(ch, :), hamming(spectral_len), spectral_overlap_len, spectral_len, fs, 'twosided')';
            Sy = cat(1, Sy, pwelch(data(ch, :), hamming(spectral_len), spectral_overlap_len, spectral_len, fs, 'twosided')');
        end
    end
end
Sy = median(Sy, 1);
Sy = Sy / sum(Sy);

params.epsilon                   = eps;
params.fs                        = 500.0;%fs;
params.filter_len                = 250;
params.lambda                    = 250.0;
params.smooth_spectrum           = true;
params.innovation_filter_type    = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE'
params.plot_results              = true;
h = equalizer_filter_design(Sx, Sy, params);


%% Test
Sx_equalized = [];
for k = 1:100%length(dataset1_fnames)
    record = fullfile(dataset1_fnames(k).folder, dataset1_fnames(k).name(1:end-4));
    record_desc = wfdb_desc(record);
    for ch = 1 : length(record_desc)
        fs = record_desc(ch).SamplingFrequency;
        data = load(record).val;
        if length(data(ch, :)) >= spectral_len && fs == 500.0
            % data(ch, :) = data(ch, :) - mean(data(ch, :));
            data(ch, :) = data(ch, :) - lp_filter_zero_phase(data(ch, :), 0.5/fs);

            data_equalized = filter(h, sqrt(sum(h.^2)), data(ch, :));
            % data_equalized = filter(h, sqrt(sum(h.^2)), randn(size(data(ch, :))));

            % Sx(k, :) = pwelch(data(ch, :), hamming(spectral_len), spectral_overlap_len, spectral_len, fs, 'twosided')';
            Sx_equalized = cat(1, Sx_equalized, pwelch(data_equalized, hamming(spectral_len), spectral_overlap_len, spectral_len, fs, 'twosided')');
        end
    end
end
Sx_equalized = mean(Sx_equalized, 1);
Sx_equalized = Sx_equalized / sum(Sx_equalized);



ff = fs * (0:spectral_len-1) / spectral_len;
figure
lgnd = {};
plot(ff - fs/2, fftshift(10*log10(abs(Sx))), 'linewidth', 3); hold on; lgnd{end+1} = 'Sx';
plot(ff - fs/2, fftshift(10*log10(abs(Sy))), 'linewidth', 3); hold on; lgnd{end+1} = 'Sy';
plot(ff - fs/2, fftshift(10*log10(abs(Sx_equalized'))), 'linewidth', 3); hold on; lgnd{end+1} = 'Sx_equalized';
grid
legend(lgnd, 'interpreter', 'none')

% time1 = (0:length(ecg1)-1)/fs;
% time2 = time1 - length(h)/2/fs;
% figure
% plot(time1, ecg1)
% hold on
% plot(time2, ecg1_equalized)
% grid
% axis tight
