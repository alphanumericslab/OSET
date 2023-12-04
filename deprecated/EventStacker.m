function [stacked_events, num_non_zeros] = EventStacker(signal, event_indexes, event_width, varargin)
% EventStacker has been deprecated. Use event_stacker instead.
warning('EventStacker has been deprecated. Use event_stacker instead.');
[stacked_events, num_non_zeros] = event_stacker(signal, event_indexes, event_width, varargin{:});