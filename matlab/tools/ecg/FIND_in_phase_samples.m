function [T0,T1] = find_in_phase_samples(peaks,phase)

L = length(peaks); % the signal length
rr_intervals = diff(find(peaks)); % the RR-intervals
% Avg, Min and Max RR-intervals
avg_rr_interval = round(mean(rr_intervals));
min_rr_interval = min(rr_intervals);
max_rr_interval = max(rr_intervals);

wlen = max_rr_interval - min_rr_interval;

T1 = zeros(1,L - avg_rr_interval - wlen);
NN = length(T1);
for t = 1 : NN
    search_window_start = max(t + avg_rr_interval - wlen, 1);
    search_window_end = min(t + avg_rr_interval + wlen,L);
   
    phase_diffs = abs(phase(t) - phase(search_window_start : search_window_end));
    [~, I] = sort(phase_diffs, 'ascend');
    for kk = 1 : length(I)
        if abs(I(kk) - t) < min_rr_interval
            continue;
        else
            T1(t) = t + avg_rr_interval + I(kk) - wlen -1;
        end
    end

end
T1 = max(T1, 1);
T1 = min(T1, NN);
T0 = 1 : NN;
