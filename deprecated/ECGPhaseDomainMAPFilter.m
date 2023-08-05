function [data_posterior_est, data_prior_est, n_var] = ECGPhaseDomainMAPFilter(data, peaks, params)
% ECGPhaseDomainMAPFilter has been deprecated. Use ecg_den_phase_domain_gp instead.
warning('ECGPhaseDomainMAPFilter has been deprecated. Use ecg_den_phase_domain_gp instead.');
[data_posterior_est, data_prior_est, n_var] = ecg_den_phase_domain_gp(data, peaks, params);