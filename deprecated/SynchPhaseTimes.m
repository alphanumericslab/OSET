function [T0,T1] = SynchPhaseTimes(peaks,phase)
% SynchPhaseTimes has been deprecated. Use synchronous_phase_samples instead.
warning('SynchPhaseTimes has been deprecated. Use synchronous_phase_samples instead.');
[T0, T1] = synchronous_phase_samples(peaks, phase);