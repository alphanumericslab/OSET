function [index, rank] = ChannelIndex14(x, ff, fs, varargin)
% ChannelIndex14 has been deprecated. Use sqi_beat_shape_variability instead.
warning('ChannelIndex14 has been deprecated. Use sqi_beat_shape_variability instead.');
if nargin == 3
    beat_width = round(0.25*fs);
else
    beat_width = varargin{1};
end
num_peak_detection_itr = 1;
ranking_mode = 'AVG_BEAT_MAX_PEAK_ALIGNED';
plot_results = 0;
[index, rank] = sqi_beat_shape_variability(x, ff, fs, method, num_peak_detection_itr, beat_width, ranking_mode, plot_results);