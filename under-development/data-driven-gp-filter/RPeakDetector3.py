# R-peak detector python implementation.
# It is based on the Pan-Tompkins R-peak detector; the "pure" P-T procedure is returing the variable "peaks". An extra step is added in order to fine tune the positions, by searching in an interval around the peaks the argmax. The "fine-tuned" positons are returned in the variable "pks".
# For testing and debugging purposes, each signal transformation is also retured (i.e. P-T involves a consecutvie series of different filters, namely Band Pass, saturation, squaring, etc)

import numpy as np
import scipy.signal as scsig
from BandPassFilter import BandPassFilter


def RPeakDetector3(signal, frac, win_l_calib, win_l_peak_search, fs):
    # print('input sig length is ' + str(len(signal)))
    k_BP = 0.7  # cut-off value
    flc = 15 / fs
    ecg_bp1 = BandPassFilter(signal, k_BP, flc)
    #
    k_BP = 0.7  # cut-off value
    fuc = 45 / fs
    ecg_bp = BandPassFilter(ecg_bp1, k_BP, fuc)
    #
    k_tan = 4
    alpha = k_tan * np.std(ecg_bp)
    ecg_tan = alpha * np.tanh(ecg_bp / alpha)
    #
    # ecg_sq = ecg_tan**2
    ecg_sq = np.abs(ecg_tan)
    #
    wind_len_MA = 0.1
    ecg_MA_sq = ecg_sq
    # ecg_MA_sq = np.sqrt(func_MA_filter(ecg_sq, wind_len_MA, fs))
    #
    maxecg = max(ecg_MA_sq)
    heightlim = frac * maxecg
    peaks, _ = scsig.find_peaks(
        ecg_MA_sq, distance=win_l_peak_search * fs, height=heightlim
    )

    if np.sum(signal[peaks] > 0) > len(signal[peaks] > 0) - np.sum(signal[peaks] > 0):
        flag = "positive"
    else:
        flag = "negative"

    #
    # print("peaks raw (as per the detector)", peaks)
    wind_len_hf = int(win_l_calib * fs)
    pks = []
    for i in range(len(peaks)):
        if i == 0:
            wind = np.arange(0, peaks[i] + wind_len_hf, 1)
            ecg_bwr_win = signal[wind]
            if flag == "positive":
                p = 0 + np.argmax(ecg_bwr_win)
            else:
                p = 0 + np.argmin(ecg_bwr_win)
            # print(i,0,peaks[i],peaks[i]+wind_len_hf,p)
        elif i == len(peaks) - 1:
            wind = np.arange(peaks[i] - wind_len_hf, len(signal), 1)
            ecg_bwr_win = signal[wind]
            if flag == "positive":
                p = peaks[i] - wind_len_hf + np.argmax(ecg_bwr_win)
            else:
                p = peaks[i] - wind_len_hf + np.argmin(ecg_bwr_win)
            # print(i,peaks[i] - wind_len_hf, peaks[i], len(signal),p)
        else:
            wind = np.arange(peaks[i] - wind_len_hf, peaks[i] + wind_len_hf, 1)
            ecg_bwr_win = signal[wind]
            if flag == "positive":
                p = peaks[i] - wind_len_hf + np.argmax(ecg_bwr_win)
            else:
                p = peaks[i] - wind_len_hf + np.argmin(ecg_bwr_win)
            # print(i,peaks[i] - wind_len_hf, peaks[i], peaks[i] + wind_len_hf,p)
        pks.append(p)
    pks = np.unique(pks)
    return (ecg_bp1, ecg_bp, ecg_tan, ecg_sq, ecg_MA_sq, peaks, np.array(pks))
