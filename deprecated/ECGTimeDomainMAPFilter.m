function [data_posterior_est, data_prior_est, n_var] = ECGTimeDomainMAPFilter(data, peaks, params)
% ECGTimeDomainMAPFilter has been deprecated. Use ecg_den_time_domain_gp instead.
warning('ECGTimeDomainMAPFilter has been deprecated. Use ecg_den_time_domain_gp instead.');
[data_posterior_est, data_prior_est, n_var] = ecg_den_time_domain_gp(data, peaks, params);