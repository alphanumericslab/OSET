function [index, rank] = sqi_level_crossing(x, th1, th2, varargin)
% sqi_level_crossing - Signal quality index based on the percentage of samples exceeding upper and lower thresholds.
%
% Usage:
%   [index, rank] = sqi_level_crossing(x, th1, th2, varargin)
%
% Inputs:
%   x: Input signal matrix with dimensions [L1, L2]
%   th1: Lower threshold vector with dimensions [L1, 1]
%   th2: Upper threshold vector with dimensions [L1, 1]
%   ranking_mode: Optional argument for sorting mode ('ascend' or 'descend')
%
% Outputs:
%   index: Signal quality index vector with dimensions [L1, 1]
%   rank: Sorted indices of the index vector based on the sorting mode
%
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version ChannelIndex0 and added sorting mode
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for sorting mode argument
if nargin > 3 && ~isempty(varargin{1})
    ranking_mode = varargin{1};
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
    I = find(x(i, :) >= th1(i) & x(i, :) <= th2(i));
    index(i) = 100 * length(I) / L2;
end

[~, rank] = sort(index, ranking_mode);