function [snr_median, snr_mean, mean_beat, median_beat, stacked_beats, noise_mean, noise_median] = events_snr(data, peak_indexes, beat_length, running_avg_len)

stacked_beats = event_stacker(data, peak_indexes, beat_length);

if nargin <= 3

    [mean_beat, ~, median_beat, ~] = robust_weighted_average(stacked_beats);

    % median_beat_replicated = ones(length(peak_indexes), 1) * median_beat;
    noise_median = stacked_beats - median_beat;%median_beat_replicated;
    snr_median = 10*log10(var(median_beat) ./ var(noise_median, [], 2));

    % mean_beat_replicated = ones(length(peak_indexes), 1) * mean_beat;
    noise_mean = stacked_beats - mean_beat;%mean_beat_replicated;
    snr_mean = 10*log10(var(mean_beat) ./ var(noise_mean, [], 2));

else
    mean_beat_moving = movmean(stacked_beats, running_avg_len, 1, "omitnan");
    mean_beat = mean(stacked_beats, 1);
    median_beat = median(stacked_beats, 1);

    noise_mean = stacked_beats - mean_beat_moving;%mean_beat_replicated;
    snr_mean = 10*log10(var(mean_beat_moving, [], 2) ./ var(noise_mean, [], 2));

    median_beat_moving = movmedian(stacked_beats, running_avg_len, 1, "omitnan");
    noise_median = stacked_beats - median_beat_moving;%mean_beat_replicated;
    snr_median = 10*log10(var(median_beat_moving, [], 2) ./ var(noise_median, [], 2));
end