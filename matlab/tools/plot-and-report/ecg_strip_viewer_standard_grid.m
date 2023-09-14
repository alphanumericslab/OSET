function ecg_strip_viewer_standard_grid(data, fs, ch_names, ref_ch, t1_small, t2_small, Title)

t1_long = 0; % beginning of the long segment in seconds
t2_long = size(data, 2)/fs; % end of long segment in seconds
num_channels = size(data, 1); % number of channels

% select multichannel grid layout
switch num_channels
    case 1
        rows = 1; cols = 1;
    case 2
        rows = 1; cols = 2;
    case 3
        rows = 1; cols = 3;
    case 4
        rows = 2; cols = 2;
    case 6
        rows = 2; cols = 3;
    case 12
        rows = 3; cols = 4;
    case 15
        rows = 3; cols = 5;
    otherwise
        cols = 4; rows = ceil(num_channels/cols);
end

% check the requested time spans
if t2_small > t2_long
    error('requested end time is longer than the signal length');
end
if t1_small < 0
    error('negative time is unacceptable');
end
    
short_ecg_lens_in_seconds = abs(t2_small - t1_small); % length of small ECGs
long_ecg_len_in_seconds = abs(t2_long - t1_long); % length of long ECG stripe

cm_per_seconds = 0.4; % equivalent width of each 1cm ECG grid square in seconds
cm_per_mV = 1.0; % equivalent height of each 1cm ECG grid square in mV
major_time_ticks_in_seconds = 0.2; % ECG vertical major time ticks in seconds
major_amp_ticks_in_mV = 0.5; % ECG horizontal amplitude major ticks in mV
% minor_time_ticks_in_seconds = 0.04; % ECG vertical minor time ticks in seconds
% minor_amp_ticks_in_mV = 0.1; % ECG horizontal amplitude minor ticks in mV
left_border_in_cm = 0.0; 
bottom_border_in_cm = 0.0;

pulse_duration = 0.6;
square_width_in_seconds = 0.2;
square_amp_in_mV = 1.0;

vertical_gap_in_mV = 4.0;

% ecg_width_in_seconds = cols * (short_ecg_lens_in_seconds + time_gap_between_plots_in_seconds) + time_gap_between_plots_in_seconds + 5.0;
min_short_ecg_total_width_in_seconds = cols * short_ecg_lens_in_seconds + (cols + 1) * major_time_ticks_in_seconds;
long_ecg_total_width_in_seconds = long_ecg_len_in_seconds + pulse_duration + 3 * major_time_ticks_in_seconds + 5.0 * square_width_in_seconds; % The 5.0 on-wards is for the pulse text
time_gap_between_plots_in_seconds = abs(long_ecg_total_width_in_seconds - cols * short_ecg_lens_in_seconds)/cols;
ecg_width_in_seconds = max(long_ecg_total_width_in_seconds, min_short_ecg_total_width_in_seconds);
side_time_gaps_long_ecg_in_seconds = abs(long_ecg_total_width_in_seconds - ecg_width_in_seconds)/2.0;

ecg_height_in_mV = (rows + 1) * vertical_gap_in_mV;
% plot_vertical_offset_in_cm = 0.0;
first_row_vertical_offset = 1.5 * vertical_gap_in_mV;
strip_row_vertical_offset = vertical_gap_in_mV/2;
all_channels_left_offset = square_width_in_seconds;

% PLOT THE FIGURES
fig = figure;
fig.Units = 'centimeters';
fig.PaperOrientation = 'Landscape';
% fig.PaperType = 'usletter';
% fig.Position = [1.0, 1.0, 26.0, 20.0];
%fig.Position = [1.0, 1.0, 40.0, 20.0];
% fig.Resize = 'off';
% fig.PaperUnits = 'centimeters';


axis_vs_fig_border_in_cm = 0.5; 
plot_width_in_cm = (ecg_width_in_seconds + 2*all_channels_left_offset)/ cm_per_seconds;
plot_height_in_cm = ecg_height_in_mV / cm_per_mV;
title_border = 1.0; % For the title place holder
fig.Position = [left_border_in_cm, bottom_border_in_cm, plot_width_in_cm, plot_height_in_cm + title_border];
% fig.PaperSize = [fig.Position(3) + 2*left_border, fig.Position(4) + 2*bottom_border];
% fig.PaperPositionMode = 'manual';

ax1 = axes;
ax1.Units = 'centimeters';
ax1.Clipping = 'on';
ax1.Position = [fig.Position(1) + axis_vs_fig_border_in_cm,...
                fig.Position(2) + axis_vs_fig_border_in_cm,...
                fig.Position(3) - 2 * axis_vs_fig_border_in_cm,...
                fig.Position(4) - 2 * axis_vs_fig_border_in_cm - title_border];
% ax1.PositionConstraint = 'innerposition';
% ax1.InnerPosition = [0, plot_vertical_offset_in_cm, plot_width_in_cm, plot_height_in_cm];
% ecg_paper_width = ax1.InnerPosition(3);
% ecg_paper_height = ax1.InnerPosition(4);


lbl = {};
t = (0 : length(data) - 1)/fs;
% mid_time_index = round(0.95 * length(t));
n1_small = find(t >= t1_small, 1, 'first');
n2_small = find(t <= t2_small, 1, 'last');
n1_long = find(t >= t1_long, 1, 'first');
n2_long = find(t <= t2_long, 1, 'last');

axes(ax1);
for ch = 1 : num_channels %size(data, 1)
    %offset = (ch - 1) * (20.0);
    %     h = subplot(rows, cols, ch);
    %     if ch == 1
    %         h_first = h;
    %     end
    %
    %     if ch == num_channels
    %         h_last = h;
    %     end
    
    row = ceil(ch / cols);
    col = ch - (row -1) * cols;
    row_reversed = rows - row + 1;
%     row = rows - floor(ch / cols)
%     col = (num_channels - ch + 1) - (row -1) * cols
    
    horizontal_offset = (col-1)*abs(short_ecg_lens_in_seconds + time_gap_between_plots_in_seconds) + all_channels_left_offset + time_gap_between_plots_in_seconds/2;
    vertical_offset = (row_reversed-1) * vertical_gap_in_mV + first_row_vertical_offset;
    
    plot(ax1, t(n1_small:n2_small) + horizontal_offset, data(ch, n1_small:n2_small) + vertical_offset, 'color', 'k', 'linewidth', 0.25); %lbl = cat(2, lbl, {'Raw signal'});
    hold on
    text(ax1, t(n1_small) + horizontal_offset, data(ch, n1_small) + vertical_offset + major_amp_ticks_in_mV, ch_names{ch}, 'color', 'b', 'fontsize', 12);
    %     hold off
    % axis tight
    % daspect(gca, [3.0, t2_small-t1_small, 1.0])
    % axis normal
end

% % % ax2 = axes;
% % % ax2.Units = 'centimeters';
% % % % ax2.Clipping = 'on';
% % % ax2.PositionConstraint = 'innerposition';
% % % ax2.InnerPosition = [0, 0.5, 20.0, 5.0]; % converting from pixels to points
% % %

% % % left_x = get(h_first, 'Position');
% % % right_x = get(h_last, 'Position');
% % %

% R-peak detection
peak_detector_params.saturate = 1;
peak_detector_params.k_sigma = 4;
peak_detector_params.hist_search_th = 0.9;
peak_detector_params.rpeak_search_wlen = 0.4; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
peak_detector_params.filter_type = 'MDMN';%'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER', 'MDMN', 'WAVELET'
[~, peak_indexes, ~] = peak_det_probabilistic(data(ref_ch, :), fs, peak_detector_params);
peak_indexes_selected_segment = peak_indexes(find(peak_indexes >= n1_long, 1, 'first') : find(peak_indexes <= n2_long, 1, 'last'));

event_width = round(1.1 * median(diff([1, peak_indexes, size(data, 2)]))); % Add the first and last indexes
if(mod(event_width, 2) == 0)
    event_width = event_width + 1;
end

stacked_beats = event_stacker(data(ref_ch, :), peak_indexes, event_width);
%         [ECG_robust_mean, ECG_robust_var, ECG_robust_median] = RWAverage(stacked_beats);
[ECG_mean, ~, ECG_median] = robust_weighted_average(stacked_beats);

lead = 0;%mean(t(peak_indexes_selected_segment(1:2)));
RR = 60.0 * fs ./ diff(peak_indexes_selected_segment);
% h01 = subplot('Position', [left_x(1), 0.21, right_x(1) - 0.18, 0.06]);
% % % plot(t(peak_indexes_selected_segment(2:end)) - lead, RR, 'k')
% % % hold on
% % % plot(t(peak_indexes_selected_segment(2:end)) - lead, RR, 'ro', 'markersize', 14)
% % % grid
% % % set(gca, 'fontsize', 14);
% % % hold off
% % % ylabel('HR (BPM)');
% % % aa = axis;
% % % aa(1) = t(n1_long) - lead;
% % % aa(2) = t(n2_long) - lead;
% % % axis(aa);
% % % % axis tight
% % % xticklabels([])

% % % h0 = subplot('Position', [left_x(1), 0.05, right_x(1) - 0.18, 0.12]);

% Plot long ECG stripe
stripe_left_offset_in_seconds = all_channels_left_offset + side_time_gaps_long_ecg_in_seconds;
plot(ax1, t(n1_long : n2_long) + stripe_left_offset_in_seconds, data(ref_ch, n1_long : n2_long) + strip_row_vertical_offset, 'color', 'k', 'linewidth', 0.25);
plot(ax1, t(peak_indexes_selected_segment) + stripe_left_offset_in_seconds, data(ref_ch, peak_indexes_selected_segment) + strip_row_vertical_offset, 'ro', 'markersize', 4)
text(ax1, t(n1_long) + stripe_left_offset_in_seconds, data(ref_ch, n1_long) + strip_row_vertical_offset + 0.5, ch_names{ref_ch}, 'color', 'b', 'fontsize', 12);

% Plot reference ECG pulse
pulse = square_amp_in_mV * (mod((0 : fs * pulse_duration - 1)/fs, 2.0*square_width_in_seconds) >= square_width_in_seconds);
pulse_time = (0:length(pulse)-1)/fs;
% plot(pulse_time + ecg_width_in_seconds - pulse_time(end), pulse + ecg_height_in_mV - max(pulse), 'k', 'linewidth', 2);
% plot(ax1, pulse_time + 15, pulse + 10.5, 'k', 'linewidth', 0.5);
% text(ax1, pulse_time(end) + 15.1, 10.5, '1mV/0.2s', 'color', 'k', 'fontsize', 6);
pulse_start_time = ceil((long_ecg_len_in_seconds + stripe_left_offset_in_seconds + major_time_ticks_in_seconds)/major_time_ticks_in_seconds)*major_time_ticks_in_seconds;
plot(ax1, pulse_time + pulse_start_time, pulse + strip_row_vertical_offset, 'k', 'linewidth', 0.5);
text(ax1, pulse_time(end) + pulse_start_time, strip_row_vertical_offset, '1mV/0.2s', 'color', 'k', 'fontsize', 12);

% % % set(ax1,'GridLineStyle','-')
set(ax1, 'fontsize', 12);
set(ax1,'LineWidth', 0.25)
% set(ax1, 'xtick', 0 : major_time_ticks_in_seconds : ecg_paper_width);
% set(ax1, 'ytick', 0 : major_amp_ticks_in_mV : ecg_paper_height);

set(ax1, 'xtick', 0 : major_time_ticks_in_seconds : (ecg_width_in_seconds + all_channels_left_offset));
set(ax1, 'ytick', 0 : major_amp_ticks_in_mV : (ecg_height_in_mV + first_row_vertical_offset));

tklen = get(ax1, 'TickLength');
set(ax1, 'TickLength', 0.1*tklen);
% % % grid minor
% % % set(ax1,'MinorGridLineStyle','--')
% % % set(ax1, 'xtick', 0 : minor_time_ticks_in_seconds : ecg_paper_width);
% % % set(ax1, 'ytick', 0 : minor_amp_ticks_in_mV : ecg_paper_height);
% % % set(ax1,'LineWidth', 0.1)
grid on
xticklabels([]);
yticklabels([]);
% axis([0, pulse_time(end) + pulse_start_time + 4*cm_per_seconds, 0, ecg_height_in_mV + square_amp_in_mV]);
axis([0, all_channels_left_offset + pulse_time(end) + pulse_start_time + 2*cm_per_seconds, 0, ecg_height_in_mV + square_amp_in_mV/2]);


if 0
    % aa = axis;
    % aa(1) = t(n1_long);
    % aa(2) = t(n2_long);
    % axis(aa);
    % axis tight
    
    avg_time_position_offset = t(n2_long) + time_gap_between_plots_in_seconds + event_width/1.5/fs;
    % % % lgnd = {};
    % % % event_width = size(stacked_beats, 2);
    tm = 1.0 * (-event_width/2 : event_width/2 - 1)/fs + avg_time_position_offset;
    % % % h1 = subplot('Position', [right_x(1), 0.05, right_x(3), 0.22]);
    plot(ax1, tm, stacked_beats' + strip_row_vertical_offset, 'color', 0.7*ones(1, 3), 'linewidth', 0.1); %lgnd = cat(2, lgnd, 'All beats');
    plot(ax1, tm, ECG_mean + strip_row_vertical_offset, 'linewidth', 0.25); %lgnd = cat(2, lgnd, 'Mean');
    plot(ax1, tm, ECG_median + strip_row_vertical_offset, 'linewidth', 0.25); %lgnd = cat(2, lgnd, 'Median');
    % % % legend(lgnd, 'Location', 'Best');
    % % % grid
    % % % chi = get(gca, 'Children'); %Returns handles to the patch and line objects
    % % % set(gca, 'Children',flipud(chi)); %Reverse the stacking order so that the averages overlay the beats
    % % % %     title(h1, 'Time-domain averaging');
    % % % xlabel('time (ms)');
    % % % % ylabel('Amplitude(mV)');
    % % % ylabel(['Lead ' ch_names{ref_ch}]);
    % % % set(gca, 'fontsize', 14)
    % % % set(gca, 'box', 'on')
    % % % aa = axis;
    % % % aa(1) = tm(1);
    % % % aa(2) = tm(end);
    % % % aa(3) = 1.1 * min([ECG_mean, ECG_median]);
    % % % aa(4) = 1.1 * max([ECG_mean, ECG_median]);
    % % % axis(aa);
end
sgtitle(Title, 'fontsize', 12)
% ax1.Units = 'centimeters';
% axis 'tight';
% aa = axis;
% % axis_borders = 1.5;
% axis(ax1, [max(0.0, aa(1) - cm_per_seconds), aa(2) + 2.0*cm_per_seconds, max(0.0, aa(3) - cm_per_mV), aa(4) + cm_per_mV])
print(fig,'LandscapePage.pdf','-dpdf', '-r0', '-bestfit')
end
