function [index, rank, yy] = sqi_power_level_crossing(x, w, th, varargin)
% Channel selection based on the number of level crossings of the signal power envelope
%
% Usage:
%   [index, rank, yy] = sqi_power_level_crossing(x, w, th, ranking_mode, normalize)
%
% Inputs:
%   x: ECG signal (rows represent channels, and columns represent samples)
%   w: window size (number of samples) for the moving average filter to compute the signal power envelope
%   th: threshold value for identifying the power envelope level crossings (% of power envelope maximum)
%   ranking_mode (optional): 'ascend' or 'descend' to sort the channels based on the calculated index. Default is 'ascend'.
%   normalize (optional): Set to 1 to return the normalized index as a percentage, or 0 to return the raw number of level crossings. Default is 1.
%
% Outputs:
%   index: a vector containing the calculated index for each channel
%   rank: a vector containing the indices that sort the channels based on their index (based on 'ranking_mode')
%   yy: matrix of power envelope signals for each channel after applying the threshold
% 
% Revision History:
%   2008: First release
%   2023: Replaced deprecated version ChannelIndex7
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check if optional arguments are provided, otherwise set default values
if nargin > 3 && ~isempty(varargin{1})
    ranking_mode = varargin{1};
    if ~isequal(ranking_mode, 'ascend') && ~isequal(ranking_mode, 'descend')
        error('Undefined ranking mode. Use ''ascend'' or ''descend''.');
    end
else
    ranking_mode = 'ascend'; % Default: sort channels in ascending order based on index
end

if nargin > 4 && ~isempty(varargin{2})
    normalize = varargin{2};
else
    normalize = 1; % Default: return normalized index as a percentage
end

L1 = size(x, 1); % Number of channels
L2 = size(x, 2); % Number of samples in each channel

yy = zeros(size(x)); % Initialize the matrix for power envelope signals after thresholding
index = zeros(L1, 1); % Initialize the vector to store the calculated index for each channel

% Calculate the power envelope for each channel and compute the index
for i = 1:L1
    s = sqrt(filter(ones(1, round(w)), round(w), x(i, :).^2)); % Calculate the power envelope using a moving average filter
    y = s - max(s) * th; % Apply the threshold to the power envelope
    yy(i, :) = y; % Store the thresholded power envelope for each channel
    
    u = y(1:end-1) .* y(2:end); % Detect zero crossings in the thresholded power envelope
    I = find(u <= 0); % Find the indices where the thresholded power envelope crosses zero

    if normalize
        index(i) = 100 * length(I) / L2; % Normalize the index as a percentage of the total number of samples
    else
        index(i) = length(I); % Return the raw number of level crossings as the index
    end
end

[~, rank] = sort(index, ranking_mode); % Sort the channels based on the calculated index
