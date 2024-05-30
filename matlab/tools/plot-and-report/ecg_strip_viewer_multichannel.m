function ecg_strip_viewer_multichannel(data, fs, varargin)
% ECG Strip Viewer - Multichannel Version
%
% A tool for visualizing multilead ECG signals.
%
% Usage:
% ecg_strip_viewer_multichannel(data, fs, ch_names, ref_ch, t1_small, t2_small, t1_long, t2_long, Title)
%
% Inputs:
%   data: Multilead ECG signal data.
%       - This is a matrix containing ECG data for multiple leads. Each row
%         corresponds to a different lead, and each column corresponds to a
%         time sample.
%
%   fs: Sampling frequency of the ECG signal.
%       - This parameter specifies the number of samples collected per second
%         for the ECG data in 'data'. It is used to convert time indices to
%         time values.
%
%   ch_names: Cell array of channel names (optional, default: empty cell array)
%       - This parameter allows you to provide names for each channel (lead) in
%         the ECG data. It should be a cell array of strings, where each string
%         corresponds to a channel name. The length of 'ch_names' should match
%         the number of rows in 'data'.
%
%   ref_ch: Reference channel for HR and median/mean plots (optional, default: 1)
%       - This parameter determines the reference channel for plotting heart rate
%         (HR) and median/mean ECG plots. It specifies the index of the reference
%         channel in the 'data' matrix.
%
%   t1_small: Start time for small plots (optional, default: 0.0)
%       - This parameter specifies the start time (in seconds) for the small
%         lead plots. It is used to define the time range for plotting the
%         individual lead signals.
%
%   t2_small: End time for small plots (optional, default: 3.0)
%       - This parameter specifies the end time (in seconds) for the small lead
%         plots. It is used to define the time range for plotting the individual
%         lead signals.
%
%   t1_long: Start time for long plots (optional, default: 0.0)
%       - This parameter specifies the start time (in seconds) for the long lead
%         plots. It is used to define the time range for plotting the reference
%         lead signal and calculating heart rate.
%
%   t2_long: End time for long plots (optional, default: 10.0)
%       - This parameter specifies the end time (in seconds) for the long lead
%         plots. It is used to define the time range for plotting the reference
%         lead signal and calculating heart rate.
%
%   Title: Title for the main figure (optional, default: empty string)
%       - This parameter allows you to provide a title for the main figure that
%         displays the multichannel ECG plots. It should be a string that
%         describes the content of the figure.
%
% Revision History:
%   2022: First release
%   2023: Renamed from deprecated version MultiLeadECGPlotter
%
% Reza Sameni, 2022-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET


if nargin > 2 && ~isempty(varargin{1})
    ch_names = varargin{1};
else
    ch_names = repmat(cell({''}), 1, size(data, 1));
end

if nargin > 3 && ~isempty(varargin{2})
    ref_ch = varargin{2};
else
    ref_ch = 1;
end

if nargin > 4 && ~isempty(varargin{3})
    t1_small = varargin{3};
else
    t1_small = 0.0;
end

if nargin > 5 && ~isempty(varargin{4})
    t2_small = varargin{4};
else
    t2_small = 3.0;
end

if nargin > 6 && ~isempty(varargin{5})
    t1_long = varargin{5};
else
    t1_long = 0.0;
end

if nargin > 7 && ~isempty(varargin{6})
    t2_long = varargin{6};
else
    t2_long = 10.0;
end

if nargin > 8 && ~isempty(varargin{7})
    Title = varargin{7};
else
    Title = '';
end

fig_left = 100;
fig_bottom = 100;
fig_width = 2000;
fig_height = 900;

% R-peak detection
% peak_detector_params.saturate = 1;
% peak_detector_params.k_sigma = 4;
% peak_detector_params.hist_search_th = 0.9;
% peak_detector_params.rpeak_search_wlen = 0.4; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
% peak_detector_params.filter_type = 'MDMN';%'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER', 'MDMN', 'WAVELET'
% [~, peak_indexes, ~] = peak_det_likelihood(data(ref_ch, :), fs, peak_detector_params);
[~, peak_indexes, ~] = peak_det_likelihood(data(ref_ch, :), fs);


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