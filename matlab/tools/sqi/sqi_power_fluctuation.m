function [index, rank, y] = sqi_power_fluctuation(x, w, varargin)
% sqi_power_fluctuation - Signal quality index based on the fluctuation of the signal power over a sliding window.
%
% Usage:
%   [index, rank, y] = sqi_power_fluctuation(x, w, ranking_mode, plot_mode)
%
% Inputs:
%   x: Input signal matrix with dimensions [L1, L2]
%   w: Sliding window length (integer value)
%   ranking_mode (optional): Sorting mode ('ascend' or 'descend'). Default: 'ascend'
%   plot_mode (optional): Plot signal and time-variant power envelope (0 or 1). Default: 0
%
% Outputs:
%   index: Signal quality index vector with dimensions [L1, 1]
%   rank: Sorted indices of the index vector based on the sorting mode
%   y: The power envelope of the input signal matrix over the sliding window of length w, with dimensions [L1, L2]
%
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version ChannelIndex2 and added sorting mode
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for sorting mode argument
if nargin > 2 && ~isempty(varargin{1})
    ranking_mode = varargin{1};
    if ~isequal(ranking_mode, 'ascend') && ~isequal(ranking_mode, 'descend')
        error('Undefined ranking mode. Please use ''ascend'' or ''descend''.');
    end
else
    ranking_mode = 'ascend';
end

% Check for plot SQI flag argument
if nargin > 3 && ~isempty(varargin{2})
    plot_sqi = varargin{2};
    if ~isequal(plot_sqi, 1) && ~isequal(plot_sqi, 0)
        error('plot_sqi flag may only be 0 or 1.');
    end
else
    plot_sqi = 0;
end

L1 = size(x, 1);

wlen = round(w);
y = zeros(size(x));
index = zeros(L1, 1);
for i = 1:L1
    channel_power = x(i, :) .^ 2;
    y(i, :) = sqrt(filter(ones(1, wlen), wlen, channel_power));
    index(i) = 100 * std(y(i, :)) / sqrt(sum(channel_power));

    if plot_sqi
        figure;
        plot(x(i, :));
        hold on;
        plot(y(i, :), 'r');
        grid;
    end
end

[~, rank] = sort(index, ranking_mode);
