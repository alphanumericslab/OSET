function [stacked_events, num_non_zeros] = EventStacker(signal, event_indexes, event_width, varargin)
%
% stacked_events = EventStacker(signal, event_indexes, event_width)
% Synchronous event stacking from input vectors and event indexes. This
% function is useful for stacking ECG beats or event related potentials
% from the R-peak locations or event times.
%
% inputs:
%   signal: the input signal in vector form
%   event_indexes: a vector or event indexes
%   event_width: the time width of the stacked events (must be odd valued)
%   method: 
%
% output:
%   stacked_events: a matrix of the form N x event_width, where N = length(event_indexes)
%   num_non_zeros: number of stacked samples per event (equal to
%   event_width except for the boundary events). Used to find the number of zero-padded samples per event
%
% Note: If event_width is even, the function updates it to the next odd value and issues a warning.
%
% Open Source Electrophysiological Toolbox, version 3.14, Nov 2020
% URL: https://gitlab.com/rsameni/OSET
%
% Copyright (C) 2020 Reza Sameni
% reza.sameni@gmail.com

if(~isvector(signal))
    error('First input should be a vector');
end

if nargin < 4
    method = 'unnormalized';
else
    method = varargin{1};
end

if(mod(event_width, 2) == 0)
    event_width = event_width + 1;
    warning('event_width must be odd valued. Automatically modified to the closest greater odd value.');
end

half_len = floor(event_width/2); % Half the window length
center_index = half_len + 1; % The center index with equal half_len samples on either side
signal_len = length(signal); % Signal length
num_events = length(event_indexes); % the number of events
stacked_events = zeros(num_events, event_width); % the matrix for stacking the events
num_non_zeros = zeros(1, num_events);

if isequal(method, 'unnormalized')
    for mm = 1 : num_events
        % The start of the event separation window (adjusted if exceeds input vector boundary)
        start = event_indexes(mm) - half_len;
        if(start < 1)
            start = 1;
        end
        left_wing_len = event_indexes(mm) - start;

        % The start of the event separation window (adjusted if exceeds input vector boundary)
        stop = event_indexes(mm) + half_len;
        if(stop > signal_len)
            stop = signal_len;
        end
        right_wing_len = stop - event_indexes(mm);

        % separate the event
        stacked_events(mm, center_index - left_wing_len : center_index + right_wing_len) = signal(start : stop);
        num_non_zeros(mm) = stop - start + 1;
    end
    % Under-development: for stretching beats (simplified dynamic time-warping)
    % elseif isequal(method, 'normalized')
    %     end_of_beat_indexes = round(event_indexes(1:end-1) + event_indexes(2:end)) / 2;
    %     for mm = 1 : length(end_of_beat_indexes)
    %         if mm > 1
    %             start = end_of_beat_indexes(mm-1) + 1;
    %         else
    %             start = 1;
    %         end
    %         stop = end_of_beat_indexes(mm);
    %         left_wing_len = event_indexes(mm) - start;
    %         right_wing_len = stop - event_indexes(mm);
    %
    %         seg_len = left_wing_len + right_wing_len + 1;
    %         % separate the event
    %         segment = zeros(1, seg_len); % the matrix for stacking the events
    %         segment(center_index - half_len : center_index + half_len) = signal(start : stop);
    %         sagment_resampled = resample(segment, event_width, seg_len);
    %         stacked_events(mm, :) = sagment_resampled;
    %         num_non_zeros(mm) = stop - start + 1;
    %     end
else
    error('undefined method');
end