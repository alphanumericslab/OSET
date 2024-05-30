function [index, rank] = sqi_beat_shape_variability(x, ff, fs, varargin)
% sqi_beat_shape_variability - Signal quality index (SQI) based on pseudo-periodicity measures
%
%   [index, rank] = sqi_beat_shape_variability(x, ff, fs, method, num_peak_det_itr, beat_width, ranking_mode, plot_results)
%
% Inputs:
%   x: Input data array (channels x samples).
%   ff: Frequency vector.
%   fs: Sampling frequency.
%   method (optional): Method for calculating the SQI. The supported methods include:
%     - 'STACKED_BEAT_EIG' (default): Calculates the SQI based on the leading eigenvalue of the stacked beat.
%     - 'AVG_BEAT_MAX': Calculates the SQI based on the maximum amplitude of the average beat.
%     - 'AVG_BEAT_MAX_PEAK_ALIGNED': Calculates the SQI based on the maximum amplitude of the average beat with R-peak amplitude alignment.
%     - 'AVG_BEAT_VAR': Calculates the SQI based on the beat variance of the average beat.
%     - 'AVG_BEAT_VAR_PEAK_ALIGNED': Calculates the SQI based on the beat variance of the average beat with R-peak amplitude alignment.
%   num_peak_det_itr (optional): Number of iterations for peak detection (default: 1).
%   beat_width (optional): Width of ECG beat used for beat stacking. Must be an odd integer (default: average heart rate time window, calculated internally).
%   ranking_mode (optional): Ranking mode for the channels ('ascend' or 'descend', default: 'descend').
%   plot_results (optional): Flag for plotting the results (default: 0).
%
% Outputs:
%   index: Signal quality index values for each channel.
%   rank: Ranked channels based on the signal quality index.
%
% Description:
%   This function computes a signal quality index (SQI) for multichannel
%   ECG signals based on different methods using the stacked-beat of the
%   ECG. The SQI is calculated by analyzing the beat shape variability in
%   the signals. The chosen method determines how the beat shape variability
%   is measured. The function returns the SQI values and the ranked channels
%   based on the SQI.
%
% Revision History:
%   2009: First release
%   2023: Added other modes and merged with deprecated versions
%       ChannelIndex9 and ChannelIndex14
%
% Reza Sameni, 2009-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if nargin > 3 && ~isempty(varargin{1})
    method = varargin{1};
else
    method = 'STACKED_BEAT_EIG';
end

if nargin > 4 && ~isempty(varargin{2})
    num_peak_det_itr = varargin{2};
else
    num_peak_det_itr = 1;
end

if nargin > 5 && ~isempty(varargin{3})
    beat_width = varargin{3};
    if mod(beat_width, 2) == 0 % beat_width must be odd
        beat_width = beat_width + 1;
        warning('beat_width must be an odd integer. updated to the next odd number.')
    end
    beat_width_auto = false;
else
    beat_width_auto = true;
end

% Check for sorting mode argument
if nargin > 6 && ~isempty(varargin{4})
    ranking_mode = varargin{4};
    if ~isequal(ranking_mode, 'ascend') && ~isequal(ranking_mode, 'descend')
        error('Undefined ranking mode. Please use ''ascend'' or ''descend''.');
    end
else
    ranking_mode = 'descend';
end

if nargin > 6 && ~isempty(varargin{4})
    plot_results = varargin{4};
else
    plot_results = 0;
end

L1 = size(x,1);
L2 = size(x,2);

index = zeros(L1,1);
mn = mean(x,2)*ones(1,L2);
x = x - mn;

for i = 1 : L1
    [peaks, peak_indexes] = peak_det_local_search(x(i,:), ff/fs, [], num_peak_det_itr);
    if beat_width_auto
        beat_width = round(mean(diff(peak_indexes)));
        if mod(beat_width, 2) == 0 % beat_width must be odd
            beat_width = beat_width + 1;
        end
    end
    switch method
        case 'STACKED_BEAT_EIG'
            stacked_beats = EventStacker(x(i, :), peak_indexes, beat_width);
            [~, D] = eig(stacked_beats*stacked_beats');
            eigs_sorted = sort(diag(D), 'descend');
            index(i) = eigs_sorted(1)/sum(eigs_sorted);
        case 'AVG_BEAT_MAX' % ChannelIndex9
            x_normalized = x(i,:) - mean(x(i,:));
            [beat_avg, beat_std] = avg_beat_calculator(x_normalized, peaks, beat_width, 'none');
            index(i) = max(abs(beat_avg)) / mean(beat_std);
        case 'AVG_BEAT_MAX_PEAK_ALIGNED' % ChannelIndex14
            x_normalized = (x(i,:) - mean(x(i,:))) / std(x(i,:));
            [beat_avg, beat_std] = avg_beat_calculator(x_normalized, peaks, beat_width, 'peak');
            index(i) = max(abs(beat_avg)) / mean(beat_std);
        case 'AVG_BEAT_VAR'
            x_normalized = x(i,:) - mean(x(i,:));
            [beat_avg, beat_std] = avg_beat_calculator(x_normalized, peaks, beat_width, 'none');
            index(i) = std(beat_avg) / mean(beat_std);
        case 'AVG_BEAT_VAR_PEAK_ALIGNED'
            x_normalized = (x(i,:) - mean(x(i,:))) / std(x(i,:));
            [beat_avg, beat_std] = avg_beat_calculator(x_normalized, peaks, beat_width, 'peak');
            index(i) = std(beat_avg) / mean(beat_std);
        otherwise
            error('Undefined method');
    end

    if plot_results
        tt = (0:length(x(i, :))-1)/fs;
        figure
        plot(tt, x(i,:));
        hold on
        plot(tt(peak_indexes), x(i, peak_indexes), 'ro', 'markersize', 14);
        grid
        xlabel('time(s)');
        ylabel('Amplitude');
        title(['channel #', num2str(i)]);
    end
end

[~, rank] = sort(index, ranking_mode);