function [baseline, baseline1] = baseline_sliding_window_twice(x, L1, L2, approach)
% baseline_sliding_window_twice - Baseline wander extraction from biomedical recordings using two stages of median or moving average filtering.
%
% Syntax: [baseline, baseline1] = baseline_sliding_window_twice(x, L1, L2, approach)
%
% Inputs:
%   x: Vector or matrix of noisy data (channels x samples).
%   L1: First stage averaging window length (in samples).
%   L2: Second stage averaging window length (in samples).
%   approach: Approach for baseline extraction.
%       - 'md': Median filtering.
%       - 'mn': Moving average.
%
% Output:
%   baseline: Vector or matrix of baseline wanders after stage 2 (channels x samples).
%   baseline1: Vector or matrix of baseline wanders after stage 1 (channels x samples).
%
%   Revision History:
%       2006: First release
%       2023: Renamed from deprecated version BaseLine1()
% 
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

N = size(x, 2);
baseline1 = zeros(size(x));
baseline = zeros(size(x));

flen1 = round(L1/2);
flen2 = round(L2/2);

if strcmp(approach, 'mn') % Moving average filter
    for j = 1:N
        index = max(j - flen1, 1) : min(j + flen1, N);
        baseline1(:, j) = mean(x(:, index), 2);
    end

    for j = 1:N
        index = max(j - flen2, 1) : min(j + flen2, N);
        baseline(:, j) = mean(baseline1(:, index), 2);
    end
elseif strcmp(approach, 'md') % Median filter
    for j = 1:N
        index = max(j - flen1, 1) : min(j + flen1, N);
        baseline1(:, j) = median(x(:, index), 2);
    end

    for j = 1:N
        index = max(j - flen2, 1) : min(j + flen2, N);
        baseline(:, j) = median(baseline1(:, index), 2);
    end
end
