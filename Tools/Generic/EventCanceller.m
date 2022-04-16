function signal_eventcancelled = EventCanceller(signal, event, event_indexes)
%
% stacked_events = EventCanceller(signal, event, event_indexes)
% Synchronous event canceller from input vectors and event indexes. This
% function is useful for stacking ECG beats or event related potentials
% from the R-peak locations or event times.
%
% inputs:
% signal: the input signal in vector form
% event: the event to be cancelled (odd lengthed)
% event_indexes: a vector or event indexes
%
% output:
% stacked_events: a matrix of the form N x event_width, where
% N = length(event_indexes)
%
% Note: If event_width is even, the function updates it to the next odd value and issues a warning.
%
% Open Source Electrophysiological Toolbox, version 3.14, Nov 2020
% URL: https://gitlab.com/rsameni/OSET
%
% Copyright (C) 2020  Reza Sameni
% reza.sameni@gmail.com

if(~isvector(signal))
    error('First input should be a vector');
end

event = event(:)';
if(mod(length(event), 2) == 0)
    event = cat(2, event, 0);
    warning('event length must be odd. Automatically appended with a 0.');
end
event_width = length(event);

half_len = floor(event_width/2); % Half the window length
center_index = half_len + 1; % The center index with equal half_len samples on either side
signal_len = length(signal); % Signal length
num_events = length(event_indexes); % the number of events
signal_eventcancelled = signal; % the signal after event cancellation

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
    
    % cancel the event
    signal_eventcancelled(start : stop) = signal(start : stop) - event(center_index - left_wing_len : center_index + right_wing_len);
end
