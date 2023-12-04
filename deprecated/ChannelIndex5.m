function [index, rank] = ChannelIndex5(x, varargin)
% ChannelIndex5 has been deprecated. Use sqi_zero_crossing instead.
warning('ChannelIndex5 has been deprecated. Use sqi_zero_crossing instead.');
[index, rank] = sqi_zero_crossing(x, varargin{:});