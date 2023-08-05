function [index, rank] = sqi_saturation_detector(x, varargin)
% sqi_saturation_detector - Signal quality index based on saturation detection.
%
% Usage:
%   [index, rank] = sqi_level_crossing(x, bins, ranking_mode, edge_bins_to_average)
%
% Inputs:
%   x: Input signal matrix with dimensions [L1, L2]
%   bins (optional): number of bins used to calculate the data histogram.
%       Default: bins = 100
%   ranking_mode (optional): Optional argument for sorting mode ('ascend'
%       or 'descend'). Default: ranking_mode = 'ascend'
%   edge_bins_to_average (optional): number of edge bins for median
%       calculation to detect saturation (excluding the first and last).
%       Default: edge_bins_to_average = 5
%
% Outputs:
%   index: Signal quality index vector with dimensions [L1, 1]
%   rank: Sorted indices of the index vector based on the sorting mode
%
% Method: Calculates the histogram of the data (per channel). If the number
%   of samples in the first or last bins are more than the median of the
%   previous few (edge_bins_to_average) edge bins, then the percentage of
%   points in the edge bins (from both tails of the distribution) are 
%   reported as the SQI. The rationale is that saturation results in more
%   than the normally expected number of samples being found in the tail of
%   the distribution.
% 
% Example 1: unsaturated
%   [index, rank] = sqi_saturation_detector(randn(10, 10000))
%
% Example 2: saturated
%   [index, rank] = sqi_saturation_detector(min(randn(10, 10000), 1))
% 
% Revision History:
%   2023: First release
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for number of bins
if nargin > 1 && ~isempty(varargin{1})
    bins = varargin{1};
else
    bins = min(100, size(x, 2));
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

% Check for the number of edge bins to average (excluding the first and
%   last)
if nargin > 3 && ~isempty(varargin{3})
    edge_bins_to_average = varargin{3};
else
    edge_bins_to_average = 5;
end


L1 = size(x, 1);
L2 = size(x, 2);

index = zeros(L1, 1);
for i = 1 : L1
    [counts, ~] = histcounts(x(i, :), bins);
    if counts(1) > median(counts(2 : 1 + edge_bins_to_average)) ...
            || counts(end) > median(counts(end - edge_bins_to_average : end - 1))
        index(i) = 100 * (counts(1) + counts(end)) / L2;
    end
end

[~, rank] = sort(index, ranking_mode);