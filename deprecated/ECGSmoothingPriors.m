function [x_filtered1, x_filtered2] = ECGSmoothingPriors(x, DiffOrder, wlen, lambda, guardlen_l, guardlen_u, adapt, withwindow)
% ECGSmoothingPriors has been deprecated. Use ecg_den_seg_wise_smoother instead.
warning('ECGSmoothingPriors has been deprecated. Use ecg_den_seg_wise_smoother instead.');
% Note: The order of outputs is different in the new version
[x_filtered2, x_filtered1] = ecg_den_seg_wise_smoother(x, DiffOrder, wlen, lambda, guardlen_l, guardlen_u, adapt, withwindow);