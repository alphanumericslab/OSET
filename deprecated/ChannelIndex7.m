function [index, rank, yy] = ChannelIndex7(x, w, th, varargin)
% ChannelIndex7 has been deprecated. Use sqi_power_level_crossing instead.
warning('ChannelIndex7 has been deprecated. Use sqi_power_level_crossing instead.');
[index, rank, yy] = sqi_power_level_crossing(x, w, th, varargin{:});
