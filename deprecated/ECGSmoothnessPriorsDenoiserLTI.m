function y = ECGSmoothnessPriorsDenoiserLTI(x, SmoothnessFactor, mode, FilterParam)
% ECGSmoothnessPriorsDenoiserLTI has been deprecated. Use ecg_den_lti_smoother instead. NOTE: ecg_den_lti_smoother accepts inputs in a different order.
warning('ECGSmoothnessPriorsDenoiserLTI has been deprecated. Use ecg_den_lti_smoother instead. NOTE: ecg_den_lti_smoother accepts inputs in a different order.');
stable_pole_mag_th = 0.9999;
y = ecg_den_lti_smoother(x, FilterParam, SmoothnessFactor, mode, stable_pole_mag_th);