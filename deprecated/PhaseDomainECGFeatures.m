function [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats_all_channels, beats_sqi] = PhaseDomainECGFeatures(data, peaks, params)
% PhaseDomainECGFeatures has been deprecated. Use ecg_features_phase_domain instead.
warning('PhaseDomainECGFeatures has been deprecated. Use ecg_features_phase_domain instead.');
[ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats_all_channels, beats_sqi] = ecg_features_phase_domain(data, peaks, params);