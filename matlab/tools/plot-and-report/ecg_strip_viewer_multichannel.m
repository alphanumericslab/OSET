function ecg_strip_viewer_multichannel(data, ch_names, fs, ref_ch, t1_small, t2_small, t1_long, t2_long, Title)
% A multilead ECG protter
%   Revision History:
%       2022: First release
%       2023: Renamed from deprecated version MultiLeadECGPlotter
%
%   Reza Sameni, 2022-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET


fig_left = 100;
fig_bottom = 100;
fig_width = 2000;
fig_height = 900;

% R-peak detection
peak_detector_params.saturate = 1;
peak_detector_params.k_sigma = 4;
peak_detector_params.hist_search_th = 0.9;
peak_detector_params.rpeak_search_wlen = 0.4; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
peak_detector_params.filter_type = 'MDMN';%'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER', 'MDMN', 'WAVELET'
[~, peak_indexes, ~] = peak_det_probabilistic(data(ref_ch, :), fs, peak_detector_params);


event_width = round(1.05 * median(diff([1, peak_indexes, size(data, 2)]))); % Add the first and last indexes
if mod(event_width, 2) == 0
    event_width = event_width + 1;
end

stacked_beats = event_stacker(data(ref_ch, :), peak_indexes, event_width);
[ECG_mean, ~, ECG_median] = robust_weighted_average(stacked_beats);


figure('Position', [fig_left, fig_bottom, fig_width, fig_height],'PaperUnits', 'points', 'PaperOrientation', 'landscape');%, 'visible','off');
%     h3 = subplot(421);
% lbl = {};
t = (0 : length(data) - 1)/fs;
% mid_time_index = round(0.95 * length(t));
n1_small = find(t >= t1_small, 1, 'first');
n2_small = find(t <= t2_small, 1, 'last');
n1_long = find(t >= t1_long, 1, 'first');
n2_long = find(t <= t2_long, 1, 'last');
peak_indexes_long = peak_indexes(find(peak_indexes >= n1_long, 1, 'first') : find(peak_indexes <= n2_long, 1, 'last'));
num_channels = size(data, 1);
switch num_channels
    case 1
        rows = 2;
        cols = 2;
    case 3
        rows = 2;
        cols = 3;
    case 6
        rows = 3;
        cols = 3;
    case 12
        rows = 4;
        cols = 4;
    case 15
        rows = 4;
        cols = 5;
end

for ch = 1 : num_channels %size(data, 1)
    %offset = (ch - 1) * (20.0);
    h = subplot(rows, cols, ch);
    if ch == 1
        h_first = h;
    end
    
    if ch == num_channels
        h_last = h;
    end
    
    plot(t(n1_small:n2_small), data(ch, n1_small:n2_small), 'color', 'k', 'linewidth', 1); %lbl = cat(2, lbl, {'Raw signal'});
    %     hold on
    %     text(t(n2_small), data(ch, n2_small) - offset - 0.3, ['lead: ', num2str(ch)], 'fontsize', 14);
    grid
    %     ylabel(['Lead ' num2str(ch)]);
    ylabel(['Lead ' ch_names{ch}]);
    %     xlabel('time (s)');
    set(gca, 'fontsize', 14);
    %     hold off
    axis tight
    % daspect(gca, [3.0, t2_small-t1_small, 1.0])
    % axis normal
end


left_x = get(h_first, 'Position');
right_x = get(h_last, 'Position');

lead = 0;%mean(t(peak_indexes_long(1:2)));
RR = 60.0 * fs ./ diff(peak_indexes_long);
subplot('Position', [left_x(1), 0.21, right_x(1) - 0.18, 0.06]);
plot(t(peak_indexes_long(2:end)) - lead, RR, 'k')
hold on
plot(t(peak_indexes_long(2:end)) - lead, RR, 'ro', 'markersize', 14)
grid
set(gca, 'fontsize', 14);
hold off
ylabel('HR (BPM)');
aa = axis;
aa(1) = t(n1_long) - lead;
aa(2) = t(n2_long) - lead;
axis(aa);
% axis tight
xticklabels([])

subplot('Position', [left_x(1), 0.05, right_x(1) - 0.18, 0.12]);
plot(t(n1_long:n2_long), data(ref_ch, n1_long:n2_long), 'color', 'k', 'linewidth', 1);
hold on
plot(t(peak_indexes_long), data(ref_ch, peak_indexes_long), 'ro', 'markersize', 14)
grid
set(gca, 'fontsize', 14);
hold off
ylabel(['Lead ' ch_names{ref_ch}]);
xlabel('time (s)');
% aa = axis;
% aa(1) = t(n1_long);
% aa(2) = t(n2_long);
% axis(aa);
axis tight

lgnd = {};
event_width = size(stacked_beats, 2);
tm = 1000.0 * (-event_width/2 : event_width/2 - 1)/fs;
subplot('Position', [right_x(1), 0.05, right_x(3), 0.22]);
plot(tm, ECG_median, 'linewidth', 3); lgnd = cat(2, lgnd, 'Median');
hold on
plot(tm, ECG_mean, 'linewidth', 3); lgnd = cat(2, lgnd, 'Mean');
plot(tm, stacked_beats', 'color', 0.7*ones(1, 3)); lgnd = cat(2, lgnd, 'All beats');
legend(lgnd, 'Location', 'Best');
grid
chi = get(gca, 'Children'); %Returns handles to the patch and line objects
set(gca, 'Children',flipud(chi)); %Reverse the stacking order so that the averages overlay the beats
%     title(h1, 'Time-domain averaging');
xlabel('time (ms)');
% ylabel('Amplitude(mV)');
ylabel(['Lead ' ch_names{ref_ch}]);
set(gca, 'fontsize', 14)
set(gca, 'box', 'on')
aa = axis;
aa(1) = tm(1);
aa(2) = tm(end);
aa(3) = 1.1 * min([ECG_mean, ECG_median]);
aa(4) = 1.1 * max([ECG_mean, ECG_median]);
axis(aa);

sgtitle(Title, 'fontsize', 18)