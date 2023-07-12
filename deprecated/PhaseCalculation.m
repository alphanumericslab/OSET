function [phase, phasepos] = PhaseCalculation(peaks)
% PhaseCalculation has been deprecated. Use phase_calculator instead.
warning('PhaseCalculation has been deprecated. Use phase_calculator instead.');
[phase, phasepos] = phase_calculator(peaks);