function [index, rank] = ChannelIndex4(x, varargin)
% ChannelIndex4 has been deprecated. Use sqi_negentropy instead.
warning('ChannelIndex4 has been deprecated. Use sqi_negentropy instead.');
[index, rank] = sqi_negentropy(x, varargin{:});