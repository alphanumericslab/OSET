function [ECGmean, ECGsd, meanPhase, ECGmedian, SamplesPerBin] = MeanECGExtraction(x,phase,bins,flag)
% MeanECGExtraction has been deprecated. Use avg_beat_calculator_phase_domain instead.
warning('MeanECGExtraction has been deprecated. Use avg_beat_calculator_phase_domain instead.');
[ECGmean, ECGsd, meanPhase, ECGmedian, SamplesPerBin] = avg_beat_calculator_phase_domain(x,phase,bins,flag);