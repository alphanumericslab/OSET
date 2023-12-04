function [stacked_events, num_non_zeros] = event_stacker(signal, event_indexes, event_bounds, varargin)
% event_stacker - Synchronous event stacking from input vectors and event indexes.
%
% Usage:
%   [stacked_events, num_non_zeros] = event_stacker(signal, event_indexes, event_bounds, method)
%
% Inputs:
%   signal: The input signal in vector form.
%   event_indexes: A vector of event indexes.
%   event_bounds:
%       - If integer scalar: event_width, the time width of the stacked events
%           (must be odd valued)
%       - if a 2-element integer vector: [left_wing_len, right_wing_len], a two
%           element vector containing the number of samples from the left and
%           right sides of each event (used for asymmetric events)
%   method (optional): Stacking method (default: 'unnormalized').
%
% Outputs:
%   stacked_events: A matrix of the form N x event_bounds, where N = length(event_indexes).
%   num_non_zeros: Number of stacked samples per event (equal to event_bounds except for the boundary events).
%                  Used to find the number of zero-padded samples per event.
%
% Note: If event_bounds is even, the function updates it to the next odd value and issues a warning.
%
%   Revision History:
%       2020: First release
%       2023: Debugged and renamed from deprecated version EventStacker
%       2023: Added asymmetric mode
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Check if the input signal is a vector
if ~isvector(signal)
    error('First input should be a vector');
end

% Determine the stacking method (default: 'unnormalized')
if nargin < 4
    method = 'unnormalized';
else
    method = varargin{1};
end

signal_len = length(signal); % Signal length
num_events = length(event_indexes); % The number of events
num_non_zeros = zeros(1, num_events);
% Separate the events using the specified method
switch method
    case 'unnormalized'
        if isscalar(event_bounds)
            % Ensure event_bounds is odd
            if mod(event_bounds, 2) == 0
                event_bounds = event_bounds + 1;
                warning('event_bounds must be odd valued (in scalar mode); automatically modified to the closest greater odd value.');
            end

            half_len = floor(event_bounds/2); % Half the window length
            center_index = half_len + 1; % The center index with equal half_len samples on either side
            stacked_events = zeros(num_events, event_bounds); % The matrix for stacking the events
            for mm = 1 : num_events
                % The start of the event separation window (adjusted if exceeds input vector boundary)
                start = event_indexes(mm) - half_len;
                if start < 1
                    start = 1;
                end
                left_wing_len = event_indexes(mm) - start;

                % The stop of the event separation window (adjusted if exceeds input vector boundary)
                stop = event_indexes(mm) + half_len;
                if stop > signal_len
                    stop = signal_len;
                end
                right_wing_len = stop - event_indexes(mm);

                % Separate the event
                stacked_events(mm, center_index - left_wing_len : center_index + right_wing_len) = signal(start : stop);
                num_non_zeros(mm) = stop - start + 1;
            end
        elseif isvector(event_bounds) && length(event_bounds) == 2
            left_wing_len0 = event_bounds(1);
            right_wing_len0 = event_bounds(2);
            intermediate_index = left_wing_len0 + 1; % The event (peak) index
            stacked_events = zeros(num_events, left_wing_len0 + right_wing_len0 + 1); % The matrix for stacking the events
            for mm = 1 : num_events
                left_wing_len = left_wing_len0;
                right_wing_len = right_wing_len0;
                % The start of the event separation window (adjusted if exceeds input vector boundary)
                start = event_indexes(mm) - left_wing_len;
                if start < 1
                    start = 1;
                end
                left_wing_len = event_indexes(mm) - start;

                % The stop of the event separation window (adjusted if exceeds input vector boundary)
                stop = event_indexes(mm) + right_wing_len;
                if stop > signal_len
                    stop = signal_len;
                end
                right_wing_len = stop - event_indexes(mm);

                % Separate the event
                stacked_events(mm, intermediate_index - left_wing_len : intermediate_index + right_wing_len) = signal(start : stop);
                num_non_zeros(mm) = stop - start + 1;
            end
        else
            error('event_bounds can be either a scalar or a two-element vector; see function help for details.')
        end

        % Under-development: for stretching beats (simplified dynamic time-warping)
        % case 'normalized'
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
        %         sagment_resampled = resample(segment, event_bounds, seg_len);
        %         stacked_events(mm, :) = sagment_resampled;
        %         num_non_zeros(mm) = stop - start + 1;
        %     end
    otherwise
        error('Undefined method');
end
