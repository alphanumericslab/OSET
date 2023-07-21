function [stacked_events, num_non_zeros] = event_stacker(signal, event_indexes, event_width, varargin)
% event_stacker - Synchronous event stacking from input vectors and event indexes.
%
% Usage:
%   [stacked_events, num_non_zeros] = event_stacker(signal, event_indexes, event_width, method)
%
% Inputs:
%   signal: The input signal in vector form.
%   event_indexes: A vector of event indexes.
%   event_width: The time width of the stacked events (must be odd valued).
%   method (optional): Stacking method (default: 'unnormalized').
%
% Outputs:
%   stacked_events: A matrix of the form N x event_width, where N = length(event_indexes).
%   num_non_zeros: Number of stacked samples per event (equal to event_width except for the boundary events).
%                  Used to find the number of zero-padded samples per event.
%
% Note: If event_width is even, the function updates it to the next odd value and issues a warning.
%
%   Revision History:
%       2020: First release
%       2023: Debugged and renamed from deprecated version EventStacker
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Check if the input signal is a vector
if(~isvector(signal))
    error('First input should be a vector');
end

% Determine the stacking method (default: 'unnormalized')
if nargin < 4
    method = 'unnormalized';
else
    method = varargin{1};
end

% Ensure event_width is odd
if(mod(event_width, 2) == 0)
    event_width = event_width + 1;
    warning('event_width must be odd valued. Automatically modified to the closest greater odd value.');
end

half_len = floor(event_width/2); % Half the window length
center_index = half_len + 1; % The center index with equal half_len samples on either side
signal_len = length(signal); % Signal length
num_events = length(event_indexes); % The number of events
stacked_events = zeros(num_events, event_width); % The matrix for stacking the events
num_non_zeros = zeros(1, num_events);

% Separate the events using the specified method
if isequal(method, 'unnormalized')
    for mm = 1 : num_events
        % The start of the event separation window (adjusted if exceeds input vector boundary)
        start = event_indexes(mm) - half_len;
        if(start < 1)
            start = 1;
        end
        left_wing_len = event_indexes(mm) - start;

        % The stop of the event separation window (adjusted if exceeds input vector boundary)
        stop = event_indexes(mm) + half_len;
        if(stop > signal_len)
            stop = signal_len;
        end
        right_wing_len = stop - event_indexes(mm);

        % Separate the event
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
    %         % Separate the event
    %         segment = zeros(1, seg_len); % The matrix for stacking the events
    %         segment(center_index - half_len : center_index + half_len) = signal(start : stop);
    %         sagment_resampled = resample(segment, event_width, seg_len);
    %         stacked_events(mm, :) = sagment_resampled;
    %         num_non_zeros(mm) = stop - start + 1;
    %     end
else
    error('Undefined method');
end
