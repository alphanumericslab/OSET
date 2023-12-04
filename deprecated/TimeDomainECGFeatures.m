function [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = TimeDomainECGFeatures(data, peaks, params)
% TimeDomainECGFeatures has been deprecated. Use ecg_features_time_domain instead.
warning('TimeDomainECGFeatures has been deprecated. Use ecg_features_time_domain instead.');
[ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = ecg_features_time_domain(data, peaks, params);