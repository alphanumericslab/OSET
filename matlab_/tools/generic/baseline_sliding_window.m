function b = baseline_sliding_window(x, L, approach)
% baseline_sliding_window - Baseline wander extraction from biomedical recordings using a single stage of median or moving average filtering.
%
% Syntax: b = baseline_sliding_window(x, L, approach)
%
% Inputs:
%   x: Vector or matrix of noisy data (channels x samples).
%   L: Averaging window length (in samples).
%   approach: Approach for baseline extraction.
%       - 'md': Median filtering.
%       - 'mn': Moving average.
%
% Output:
%   b: Vector or matrix of baseline wanders (channels x samples).
%
%   Revision History:
%       2006: First release
%       2023: Renamed from deprecated version BaseLine1()
% 
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET


N = size(x, 2);
b = zeros(size(x));
flen = floor(L/2);

if strcmp(approach, 'mn')
    % Moving average filter
    for j = 1:N
        index = max(j - flen, 1) : min(j + flen, N);
        b(:, j) = mean(x(:, index), 2);
    end
elseif strcmp(approach, 'md')
    % Median filter
    for j = 1:N
        index = max(j - flen, 1) : min(j + flen, N);
        b(:, j) = median(x(:, index), 2);
    end
end
