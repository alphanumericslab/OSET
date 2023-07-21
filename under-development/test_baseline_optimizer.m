% test baseline wander removal by periodicity maximization
% Status: UNDER TEST
% Reza Sameni, 2023

clear;
close all;

% load patient165_s0323lre; data = data(1:25000, 2)'; fs = 1000;
load SampleECG128Hz; data = y(:, 3)' ; fs = 128;

% data = data - lp_filter_zero_phase(data, 0.05/fs);

data_hp = data - lp_filter_zero_phase(data, 2.0/fs);
num_rounds = 3;
hr_update_fraction = 1.05;
omit_close_peaks = true;
[peaks, peak_indexes] = peak_detection_local_search(data_hp, 1.0/fs, [], num_rounds, hr_update_fraction, omit_close_peaks);


mn_wlens = 0.2 : 0.05 : 2.0;
md_wlens = 0.2 : 0.05 : 1.5;
sqi_index = zeros(length(md_wlens), length(mn_wlens));
for i = 1 : length(md_wlens)
    for j = 1 : length(mn_wlens)
        % baseline = baseline_sliding_window(data, round(fs * mn_wlens(k)), method);
        % baseline = baseline_sliding_window(data, round(fs * 1.0), 'md');

        %         bl = baseline_sliding_window(data, round(fs * md_wlens(i)), 'md');
        bl = baseline_sliding_window_twice(data, round(fs * md_wlens(i)), round(fs * md_wlens(i)), 'md');
        baseline = baseline_sliding_window(bl, round(fs * mn_wlens(j)), 'mn');

        ecg = data - baseline;

        event_width = round(mean(diff(peak_indexes)));
        if mod(event_width, 2) == 0
            event_width = event_width + 1;
        end
        ecg_stacked = EventStacker(ecg, peak_indexes, event_width);
        baseline_stacked = EventStacker(baseline, peak_indexes, event_width);
        % baseline_stacked = EventStacker(data, peak_indexes, event_width);

        ecg_beat_median = median(ecg_stacked, 1);
        %         ecg_stacked_demeaned = ecg_stacked - ecg_beat_median;
        %         ecg_beat_std = std(ecg_stacked_demeaned, [], 1);
        ecg_beat_std = std(ecg_stacked, [], 1);
        ecg_snr = mean(ecg_beat_median.^2) / mean(ecg_beat_std.^2);

        baseline_beat_mean = mean(baseline_stacked, 1);
        %         baseline_stacked_demeaned = baseline_stacked - baseline_beat_mean;
        %         baseline_beat_std = std(baseline_stacked_demeaned, [], 1);
        baseline_beat_std = std(baseline_stacked, [], 1);
        baseline_snr = mean(baseline_beat_mean.^2) / mean(baseline_beat_std.^2);

        baseline_most_periodic_peak = max(mean(baseline_stacked.^2, 1)./var(baseline_stacked, [], 1));

        %         cov_ecg = ecg_stacked * ecg_stacked';
        %         cov_baseline = baseline_stacked * baseline_stacked';
        %                         cov_ecg = ecg_stacked' * ecg_stacked;
        %                         cov_baseline = baseline_stacked' * baseline_stacked;
        %         sqi_index(i, j) = eigs(cov_ecg, cov_baseline, 1);
        %         sqi_index(i, j) = trace(cov_ecg)/trace(cov_baseline);

        %         sqi_index(i, j) = ecg_snr / baseline_snr;
        %         sqi_index(i, j) = ecg_snr / baseline_most_periodic_peak;
        %         sqi_index(i, j) = ecg_snr / min(baseline_beat_std);
        %         sqi_index(i, j) = mean(ecg_beat_median.^2) / mean(baseline_beat_mean.^2);

        %         baseline_stacked_lf = zeros(size(baseline_stacked));
        %         for mm = 1 : size(baseline_stacked, 1)
        %             baseline_stacked_lf(mm, :) = wden(baseline_stacked(mm, :),'rigrsure','s','sln', 1, 'coif5');
        %         end
        %         baseline_stacked_lf = baseline_stacked_lf - mean(baseline_stacked_lf, 2);
        baseline_stacked_lf = baseline_stacked;
        baseline_stacked_lf_median = median(baseline_stacked_lf, 1);

        %         sqi_index(i, j) = (size(ecg_stacked, 1) * mean(ecg_beat_median.^2) / trace(ecg_stacked*ecg_stacked') ) / ...
        %             (size(baseline_stacked_lf, 1) * mean(baseline_stacked_lf_median.^2) / trace(baseline_stacked_lf*baseline_stacked_lf') );

        %         sqi_index(i, j) = (size(ecg_stacked, 1) * mean(ecg_beat_median.^2) / trace(ecg_stacked*ecg_stacked') ) / ...
        %             (size(baseline_stacked_lf, 1) * mean(baseline_stacked_lf_median.^2) / trace(cov(baseline_stacked_lf')) );

        sqi_index(i, j) = (size(ecg_stacked, 1) * mean(ecg_beat_median.^2) / trace(ecg_stacked*ecg_stacked') ) / ...
                mean(baseline_stacked_lf_median.^2);
        %         sqi_index(i, j) = size(ecg_stacked, 1) * mean(ecg_beat_median.^2) / trace(ecg_stacked*ecg_stacked');

        disp([i, j])
    end
end

[~, opt_index] = max(sqi_index(:));
[opt_md_index, opt_mn_index] = ind2sub(size(sqi_index), opt_index);
bl = baseline_sliding_window(data, round(fs * md_wlens(opt_md_index)), 'md');
baseline_opt = baseline_sliding_window(bl, round(fs * mn_wlens(opt_mn_index)), 'mn');
ecg_opt = data - baseline_opt;

ecg_opt_stacked = EventStacker(ecg_opt, peak_indexes, event_width);
baseline_opt_stacked = EventStacker(baseline_opt, peak_indexes, event_width);

figure
subplot(121)
plot(ecg_opt_stacked');
hold on
plot(std(ecg_opt_stacked, [], 1), 'k', 'linewidth', 2);
grid
title('ecg_opt_stacked', 'interpreter', 'none');
subplot(122)
plot(baseline_opt_stacked');
hold on
plot(std(baseline_opt_stacked, [], 1), 'k', 'linewidth', 2);
grid
title('baseline_opt_stacked', 'interpreter', 'none');


figure
plot(mn_wlens, sqi_index);
grid
title('sqi_index', 'interpreter', 'none');
xlabel('mn_wlens', 'interpreter', 'none');

t = (0 : length(data)-1)/fs;
lgnd = {};
figure
plot(t, data); lgnd = cat(1, lgnd, 'raw signal');
hold on
plot(t(peak_indexes), data(peak_indexes), 'ro', 'markersize', 18); lgnd = cat(1, lgnd, 'R-peaks');
plot(t, baseline_opt); lgnd = cat(1, lgnd, 'baseline opt');
plot(t, ecg_opt); lgnd = cat(1, lgnd, 'baseline removed');
grid
legend(lgnd)

figure
[md_wlens_mesh, mn_wlens_mesh] = meshgrid(md_wlens, mn_wlens);
surf(md_wlens_mesh, mn_wlens_mesh, sqi_index');
xlabel('md_wlens_mesh', 'interpreter', 'none');
ylabel('mn_wlens_mesh', 'interpreter', 'none');
zlabel('sqi_index', 'interpreter', 'none');

