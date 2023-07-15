function [score1, score2, peaks1, peaks2, rank1, rank2] = ChannelIndex12(x,f,fs, varargin)
% ChannelIndex12 has been deprecated. Use sqi_multi_matched_filter instead.
warning('ChannelIndex12 has been deprecated. Use sqi_multi_matched_filter instead.');
num_peak_detection_itr = 1;
template_type = 0;
template_params = [];
plot_templates = 0;
[score1, score2, peaks1, peaks2, rank1, rank2] = sqi_multi_matched_filter(x, f, fs, num_peak_detection_itr, template_type, template_params, plot_templates);