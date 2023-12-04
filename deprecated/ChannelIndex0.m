function [index, rank] = ChannelIndex0(x, th1, th2, varargin)
% ChannelIndex0 has been deprecated. Use sqi_level_crossing instead.
warning('ChannelIndex0 has been deprecated. Use sqi_level_crossing instead.');
[index, rank] = sqi_level_crossing(x, th1, th2, varargin{:});