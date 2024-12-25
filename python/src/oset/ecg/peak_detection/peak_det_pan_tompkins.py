import argparse
import math
import numpy as np
from scipy.signal import butter, filtfilt, convolve, resample, lfilter


def peak_det_pan_tompkins(
    ecg_data, fs, ecg_polarity=None, qrs_width=0.150, refracT=0.360
):
    """
    peak_det_pan_tompkins - R-peak detector based on Pan-Tompkins method.
    QRS Detection with Pan-Tompkins algorithm.
    Pan et al. 1985: A Real-Time QRS Detection Algorithm.

    This function implements the Pan-Tompkins algorithm for R-peak detection in ECG signals,
    with a simplified post-detection R-peak selection logic

    Args:
        ecg_data (numpy.ndarray): Vector of input ECG data.
        fs (int): Sampling rate in Hz.
        ecg_polarity (numpy.ndarray, optional): Search for positive (flag=1) or negative (flag=0) peaks.
            By default, the skewness value of the bandpass-filtered signal determines the peak sign.
        qrs_width (float, optional): Expected maximal length of QRS complexes in seconds.
            Default is 0.150 seconds (150 ms).
        refracT (float, optional): Refractory time for T-wave in seconds.
            Default is 0.360 seconds (360 ms).

    Reference:
        Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
        Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532

    Revision History:
        July 2023: Translated to Python from Matlab (peak_det_pan_tompkins.m)
        June 2024: Updated to match the new MATLAB version (peak_det_pan_tompkins.m)

    Amulya Jain, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    raise NotImplementedError(
        "This function has not been implemented yet. It is a work in progress."
    )
    # Default parameters
    if ecg_polarity is None:
        ecg_polarity = np.array([])

    data = np.copy(ecg_data)

    fs_pt = 200
    total_delay = 0

    g = math.gcd(fs, fs_pt)
    params = {}

    num_samples = len(data) * int(fs_pt / g) // int(fs / g)

    data = resample(data, num_samples)

    # Pan-Tomkins signal processing

    # 1 remove the mean
    data = data - np.mean(data)

    # 2 low-pass filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
    b = np.array([1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1])
    a = np.array([1, -2, 1])
    filtered_data = lfilter(b, a, data)
    total_delay += 6

    # 3 high-pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
    b = np.zeros(33)
    b[0] = -1
    b[16] = 32
    b[32] = 1
    a = np.array([1, 1])
    filtered_data = lfilter(b, a, filtered_data)
    total_delay += 16
    filter_delay = total_delay

    # Differentiation to enhance R-peaks
    # H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))

    b = np.array([1, 2, 0, -2, -1]) / 8 * fs

    diff_data = lfilter(b, 1, filtered_data)
    total_delay += 2

    # Squaring to further emphasize R-peaks
    squared_data = diff_data**2

    # Moving average integration
    window_length = np.round(qrs_width * fs)

    integrated_data = np.convolve(
        squared_data, np.array([window_length, window_length])
    )

    # remove delays & resample
    integrated_data = np.concatenate(
        (
            integrated_data[total_delay:],
            np.repeat(integrated_data[-1], total_delay - 1, axis=0),
        ),
        axis=0,
    ).T
    integrated_data = resample(integrated_data, int(fs / g), int(fs_pt / g))
    integrated_data = integrated_data / np.std(integrated_data)

    filtered_data = np.concatenate(
        (
            filtered_data[filter_delay:],
            np.repeat(filtered_data[-1], filter_delay - 1, axis=0),
        ),
        axis=0,
    ).T
    filtered_data = resample(filtered_data, int(fs / g), int(fs_pt / g))

    if len(ecg_polarity) == 0:
        filtered_data = filtered_data * np.sign(skew(filtered_data))
    else:
        filtered_data = filtered_data * np.sign(ecg_polarity - 0.5)

    # High-pass filter
    b_hp, a_hp = butter(5, params["fc_low"] / (fs / 2), "high")
    filtered_data_hp = filtfilt(
        b_hp, a_hp, ecg_data, padlen=3 * (max(len(b_hp), len(a_hp)) - 1), axis=-1
    )

    # Low-pass filter
    b_lp, a_lp = butter(5, params["fc_high"] / (fs / 2), "low")
    filtered_data = filtfilt(
        b_lp,
        a_lp,
        filtered_data_hp,
        padlen=3 * (max(len(b_lp), len(a_lp)) - 1),
        axis=-1,
    )

    # Differentiation to enhance R-peaks
    diff_data = np.diff(filtered_data, prepend=0)

    # Squaring to further emphasize R-peaks
    squared_data = diff_data**2

    # Moving average integration
    window_length = round(params["window_length"] * fs)
    window = np.ones(window_length) / window_length
    integrated_data = convolve(squared_data, window, mode="same")

    # Amplitude thresholding
    threshold = params["threshold_ratio"] * np.max(integrated_data)

    # Refractory period to avoid detecting multiple peaks within a short duration
    refractory_half_period = round(params["refractory_period"] * fs / 2)

    # Find the local peaks that satisfy both amplitude and width condition
    peaks = np.zeros(ecg_data.shape[0])
    for k in range(1, ecg_data.shape[0] + 1):
        search_win_start = max(1, k - refractory_half_period)
        search_win_stop = min(ecg_data.shape[0], k + refractory_half_period + 1)
        if (
            threshold
            < integrated_data[k - 1]
            == np.max(integrated_data[search_win_start:search_win_stop])
        ):
            peaks[k - 2] = 1
    # Get peak indexes
    peak_indexes = np.nonzero(peaks)[0] + 1

    return peaks, peak_indexes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_det_pan_tompkins - R-peak detector based on Pan-Tompkins method.
    QRS Detection with Pan-Tompkins algorithm.
    Pan et al. 1985: A Real-Time QRS Detection Algorithm.

    This function implements the Pan-Tompkins algorithm for R-peak detection in ECG signals,
    with a simplified post-detection R-peak selection logic

    Args:
        ecg_data (numpy.ndarray): Vector of input ECG data.
        fs (int): Sampling rate in Hz.
        ecg_polarity (numpy.ndarray, optional): Search for positive (flag=1) or negative (flag=0) peaks.
            By default, the skewness value of the bandpass-filtered signal determines the peak sign.
        qrs_width (float, optional): Expected maximal length of QRS complexes in seconds.
            Default is 0.150 seconds (150 ms).
        refracT (float, optional): Refractory time for T-wave in seconds.
            Default is 0.360 seconds (360 ms).

    Reference:
        Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
        Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532

    Revision History:
        July 2023: Translated to Python from Matlab (peak_det_pan_tompkins.m)
        June 2024: Updated to match the new MATLAB version (peak_det_pan_tompkins.m)

    Amulya Jain, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
