function [index, rank] = sqi_negentropy(x, varargin)
% sqi_negentropy - Signal quality index based on negentropy estimation.
%
% Usage:
%   [index, rank] = sqi_negentropy(x, ranking_mode)
%
% Inputs:
%   x: Input signal matrix with dimensions [L1, L2]
%   ranking_mode: Optional argument for sorting mode ('ascend' or 'descend'). Default: 'ascend'
%
% Outputs:
%   index: Channel index vector with dimensions [L1, 1]
%   rank: Sorted indices of the index vector based on the sorting mode
%
% Reference:
%   A. Hyvärinen, J. Karhunen, and E. Oja, Independent Component Analysis.
%   Wiley-Interscience, 2001, Chapter 5
%
% Revision History:
%   2008: First release
%   2023: Corrected fnegentropy definition according to the reference and renamed from deprecated version ChannelIndex4 and added sorting mode
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for sorting mode argument
if nargin > 1 && ~isempty(varargin{1})
    ranking_mode = varargin{1};
    if ~isequal(ranking_mode, 'ascend') && ~isequal(ranking_mode, 'descend')
        error('Undefined ranking mode. Please use ''ascend'' or ''descend''.');
    end
else
    ranking_mode = 'ascend';
end

L1 = size(x, 1);

index = zeros(L1, 1);

for i = 1:L1
    % Channel centralization & normalization
    x(i, :) = (x(i, :) - mean(x(i, :))) / std(x(i, :));

    switch method
        case 'J' % See eq. (5.35) in Hyvärinen et al. 2001
            index(i) = mean(x(i, :).^3)^2 / 12 + kurtosis(x(i, :))^2 / 48;
        case 'Ja' % See eq. (5.47) in Hyvärinen et al. 2001
            k1 = 36 / (8 * sqrt(3) - 9);
            k2a = 1 / (2 - 6 / pi);
            index(i) = k1 * mean(x(i, :) .* exp(-x(i, :).^2 / 2))^2 + k2a * (mean(abs(x(i, :))) - sqrt(2 / pi))^2;
        case 'Jb' % See eq. (5.48) in Hyvärinen et al. 2001
            k1 = 36 / (8 * sqrt(3) - 9);
            k2b = 24 / (16 * sqrt(3) - 27);
            index(i) = k1 * mean(x(i, :) .* exp(-x(i, :).^2 / 2))^2 + k2b * (mean(exp(-x(i, :).^2 / 2)) - sqrt(1 / 2))^2;
    end
end

[~, rank] = sort(index, ranking_mode);
