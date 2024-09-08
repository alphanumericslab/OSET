import numpy as np
from scipy.signal import medfilt
from scipy.stats import zscore

def em_interval_calc(ecg_rpeaks_index, thr_max=0.5, thr_fast_max=0.4, max_lag=20):
#    if ecg_rpeaks_index.shape[1] == 1:
#        rr_intervals_ecg = np.diff(ecg_rpeaks_index.flatten())
#    elif ecg_rpeaks_index.shape[1] == 2:
#        rr_intervals_ecg = ecg_rpeaks_index[:, 1] - ecg_rpeaks_index[:, 0]
    rr_intervals_ecg = np.diff(ecg_rpeaks_index)

    md_rr_intervals_ecg = medfilt(rr_intervals_ecg, kernel_size=2*max_lag+1)
    diff_outliers = rr_intervals_ecg / md_rr_intervals_ecg - 1
    ind_nan_rr = np.abs(diff_outliers) > thr_max
    rr_intervals_ecg[ind_nan_rr] = md_rr_intervals_ecg[ind_nan_rr]

    if np.any(ind_nan_rr):
        md_rr_intervals_ecg = medfilt(rr_intervals_ecg, kernel_size=2*int(np.ceil(max_lag/3))+1)
        diff_outliers = rr_intervals_ecg / md_rr_intervals_ecg - 1
        ind_nan_rr = np.abs(diff_outliers) > thr_fast_max
        rr_intervals_ecg[ind_nan_rr] = md_rr_intervals_ecg[ind_nan_rr]

        diff_outliers = rr_intervals_ecg - md_rr_intervals_ecg
        std_rr_intervals_ecg = np.array([np.std(diff_outliers[max(0, i-100):min(len(diff_outliers), i+101)]) for i in range(len(diff_outliers))])
        z_score_rr = np.abs(diff_outliers / std_rr_intervals_ecg)
        ind_nan_zrr = z_score_rr > 5
        rr_intervals_ecg[ind_nan_zrr] = md_rr_intervals_ecg[ind_nan_zrr]

    # Linear interpolation for missing values
    nan_mask = np.isnan(rr_intervals_ecg)
    rr_intervals_ecg[nan_mask] = np.interp(np.flatnonzero(nan_mask), np.flatnonzero(~nan_mask), rr_intervals_ecg[~nan_mask])

    return rr_intervals_ecg