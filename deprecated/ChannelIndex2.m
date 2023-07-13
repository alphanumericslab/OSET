function [index, rank, y] = ChannelIndex2(x, w, varargin)
% ChannelIndex2 has been deprecated. Use sqi_power_fluctuation instead.
warning('ChannelIndex2 has been deprecated. Use sqi_power_fluctuation instead.');
[index, rank, y] = sqi_power_fluctuation(x, w, varargin{:});