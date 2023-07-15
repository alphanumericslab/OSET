function [index, rank] = sqi_pseudo_periodicity(x,ff,fs, varargin)
% Signal quality index based on pseudo-periodicty measures
%
% Revision History:
%   2019: First release
%   2023: Added other modes and renamed from deprecated version ChannelIndex13
%
% References:
%     - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2009). A deflation procedure for subspace decomposition. IEEE Transactions on Signal Processing, 58(4), 2363-2374.
%     - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2008). Multichannel electrocardiogram decomposition using periodic component analysis. IEEE transactions on biomedical engineering, 55(8), 1935-1940.
%
% Reza Sameni, 2019-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if nargin > 3 && ~isempty(varargin{1})
    method = varargin{1};
else
    method = 'EVD';
end

if nargin > 4 && ~isempty(varargin{2})
    num_peak_detection_itr = varargin{2};
else
    num_peak_detection_itr = 1;
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
L2 = size(x,2);

index = zeros(L1,1);
mn = mean(x,2)*ones(1,L2);
x = x - mn;

for i = 1:L1
    [~, peak_indexes] = peak_detection_local_search(x(i,:), ff/fs, [], num_peak_detection_itr);
    event_width = round(mean(peak_indexes));
    if mod(event_width, 2) == 0
        event_width = event_width + 1;
    end
    stacked_beats = EventStacker(x(i, :), peak_indexes, event_width);
    [~, D] = eig(stacked_beats*stacked_beats');
    eigs_sorted = sort(diag(D), 'descend');
    index(i) = eigs_sorted(1)/sum(eigs_sorted);

%{
    [tt0, tt1] = SynchPhaseTimes2(peaks);
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
            [~, D] = eig(A,B);
            eigs_sorted = sort(diag(D), 'descend');
            index(i) = eigs_sorted(1)/sum(eigs_sorted);
        case 'TRACE'
            index(i) = trace(A)/trace(B);
        case 'STACKED_BEATS'
            event_width = mean(peak_indexes);
            if mod(event_width, 2) == 0
                event_width = event_width + 1;
            end
            stacked_beats = EventStacker(x(i, :), peak_indexes, event_width);
            [~, D] = eig(stacked_beats*stacked_beats');
            eigs_sorted = sort(diag(D), 'descend');
            index(i) = eigs_sorted(1)/sum(eigs_sorted);
        otherwise
            error('undefined method');
    end
%}

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
end

[~, rank] = sort(index, ranking_mode);
