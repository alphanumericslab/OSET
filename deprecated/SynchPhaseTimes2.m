function [T0, T1] = SynchPhaseTimes2(peaks)
% SynchPhaseTimes2 has been deprecated. Use synchronous_phase_samples instead.
warning('SynchPhaseTimes2 has been deprecated. Use synchronous_phase_samples instead.');
[T0, T1] = synchronous_phase_samples(peaks);