function [index, rank] = sqi_pseudo_periodicity(x, ff, fs, varargin)
% Signal quality index based on pseudo-periodicity measures
%
% [index, rank] = sqi_pseudo_periodicity(x, ff, fs, method, num_peak_det_itr, ranking_mode, plot_results)
%
% This function calculates the signal quality index (SQI) based on
%   pseudo-periodicity measures of the input data.
%
% Inputs:
%   x: Matrix of input data (channels x samples)
%   ff: Fundamental frequency of the signal in Hz
%   fs: Sampling rate of the signal in Hz
%   method (optional): Method for SQI calculation ('STACKED_BEATS', 'EVD', 'GEVD', 'TRACE'). Default is 'STACKED_BEATS'
%   num_peak_det_itr (optional): Number of iterations for peak detection (optional, default is 1)
%   ranking_mode (optional): Sorting mode for the output ranking ('ascend' or 'descend', optional, default is 'descend')
%   plot_results (optional): Flag to plot the results (0 or 1, optional, default is 0)
%
% Outputs:
%   index: Vector of SQI values for each channel
%   rank: Indices that sort the channels based on SQI values
%
% Revision History:
%   2019: First release
%   2023: Added other modes and renamed from deprecated version ChannelIndex13
% 
% References:
%   - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2009). A deflation
%       procedure for subspace decomposition. IEEE Transactions on Signal
%       Processing, 58(4), 2363-2374.
%   - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2008). Multichannel
%       electrocardiogram decomposition using periodic component analysis. IEEE
%       transactions on biomedical engineering, 55(8), 1935-1940.
%
% Reza Sameni, 2019-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for optional input arguments and set default values
if nargin > 3 && ~isempty(varargin{1})
    method = varargin{1};
else
    method = 'STACKED_BEATS';
end

if nargin > 4 && ~isempty(varargin{2})
    num_peak_det_itr = varargin{2};
else
    num_peak_det_itr = 1;
end

% Check for sorting mode argument
if nargin > 5 && ~isempty(varargin{3})
    ranking_mode = varargin{3};
    if ~isequal(ranking_mode, 'ascend') && ~isequal(ranking_mode, 'descend')
        error('Undefined ranking mode. Please use ''ascend'' or ''descend''.');
    end
else
    ranking_mode = 'descend';
end

if nargin > 6 && ~isempty(varargin{4})
    plot_results = varargin{4};
else
    plot_results = 0;
end

L1 = size(x,1);

index = zeros(L1, 1);
x = x - mean(x,2);

for i = 1 : L1
    [peaks, peak_indexes] = peak_det_local_search(x(i,:), ff/fs, [], num_peak_det_itr);
    if plot_results
        tt = (0:length(x(i, :))-1)/fs;
        figure
        plot(tt, x(i,:));
        hold on
        plot(tt(peak_indexes), x(i, peak_indexes), 'ro', 'markersize', 14);
        grid
        xlabel('time(s)');
        ylabel('Amplitude');
        title(['channel #', num2str(i)]);
    end

    [tt0, tt1] = synchronous_phase_samples(peaks);
    % periodic component analysis stage
    A = x(:,tt0)*x(:,tt1)'/length(tt0);
    B = x(:,tt0)*x(:,tt0)'/length(tt0);

    A = (A + A') / 2;
    B = (B + B') / 2;

    switch method
        case 'EVD'
            [~, D] = eig(B);
            eigs_sorted = sort(diag(D), 'descend');
            index(i) = eigs_sorted(1)/sum(eigs_sorted);
        case 'GEVD'
            [~, D] = eig(A, B);
            eigs_sorted = sort(diag(D), 'descend');
            index(i) = eigs_sorted(1)/sum(eigs_sorted);
        case 'TRACE'
            index(i) = trace(A)/trace(B);
        case 'STACKED_BEATS'
            event_width = mean(peak_indexes);
            if mod(event_width, 2) == 0
                event_width = event_width + 1;
            end
            stacked_beats = event_stacker(x(i, :), peak_indexes, event_width);
            [~, D] = eig(stacked_beats*stacked_beats');
            eigs_sorted = sort(diag(D), 'descend');
            index(i) = eigs_sorted(1)/sum(eigs_sorted);
        otherwise
            error('undefined method');
    end

end

[~, rank] = sort(index, ranking_mode);
end
