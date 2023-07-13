function [index, rank] = sqi_zero_crossing(x, varargin)
% sqi_zero_crossing - Signal quality index based on the percentage of zero-crossings.
%
% Usage:
%   [index, rank] = sqi_zero_crossing(x, remove_mean)
%
% Inputs:
%   x: Input signal matrix with dimensions [L1, L2]
%   remove_mean: Optional arguments for remove_mean flag (0 or 1) and sorting mode ('ascend' or 'descend')
%
% Outputs:
%   index: Signal quality index vector with dimensions [L1, 1]
%   rank: Sorted indices of the index vector based on the sorting mode
%
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version ChannelIndex5 and added sorting mode
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for remove_mean flag argument
if nargin > 1 && ~isempty(varargin{1})
    remove_mean = varargin{1};
    if ~isequal(remove_mean, 0) && ~isequal(remove_mean, 1)
        error('''remove_mean'' should be 0 or 1');
    end
else
    remove_mean = 1;
end

% Check for sorting mode argument
if nargin > 2 && ~isempty(varargin{2})
    ranking_mode = varargin{2};
    if ~isequal(ranking_mode, 'ascend') && ~isequal(ranking_mode, 'descend')
        error('Undefined ranking mode. Please use ''ascend'' or ''descend''.');
    end
else
    ranking_mode = 'ascend';
end

L1 = size(x, 1);
L2 = size(x, 2);

index = zeros(L1, 1);

for i = 1:L1
    if remove_mean
        x(i, :) = x(i, :) - mean(x(i, :));
    end
    sgn = x(i, 1:end-1) .* x(i, 2:end);
    I = find(sgn <= 0);
    index(i) = 100 * length(I) / (L2 - 1);
end

[~, rank] = sort(index, ranking_mode);
