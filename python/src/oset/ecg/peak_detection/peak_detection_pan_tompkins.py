import argparse

import numpy as np
from scipy.signal import butter, filtfilt, convolve


def peak_detection_pan_tompkins(data, fs, *args):
    """
    peak_detection_pan_tompkins - R-peak detector based on Pan-Tompkins method.

    Args:
        data (numpy.ndarray): Vector of input ECG data
        fs (float): Sampling rate in Hz
        fc_low (float, optional): BP filter lower cutoff frequency in Hz (default: 5.0 Hz)
        fc_high (float, optional): BP filter upper cutoff frequency in Hz (default: 15.0 Hz)
        window_length (float, optional): Integration window length in seconds (default: 0.150 s)
        threshold_ratio (float, optional): Threshold ratio for peak detection (default: 0.2)
        refractory_period (float, optional): Refractory period in seconds (default: 0.2 s)

    Returns:
        peaks (numpy.ndarray): Vector of R-peak impulse train
        peak_indexes (numpy.ndarray): Vector of R-peak indexes

    References:
        Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
        Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532

    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_pan_tompkins.m)

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    default_params = {
        'fc_low': 5.0,
        'fc_high': 15,
        'window_length': 0.150,
        'threshold_ratio': 0.2,
        'refractory_period': 0.200
    }

    def fill_default_params(params, default_params):
        """Helper function to fill missing fields in params struct with default values"""
        for key, value in default_params.items():
            if key not in params:
                params[key] = value
        return params

    # Parse input arguments
    params = default_params.copy()
    if len(args) > 0 and isinstance(args[0], dict):
        # params dict is provided
        params = fill_default_params(args[0], default_params)

    # High-pass filter
    b_hp, a_hp = butter(5, params['fc_low'] / (fs / 2), 'high')
    filtered_data_hp = filtfilt(b_hp, a_hp, data, padlen=3 * (max(len(b_hp), len(a_hp)) - 1), axis=-1)

    # Low-pass filter
    b_lp, a_lp = butter(5, params['fc_high'] / (fs / 2), 'low')
    filtered_data = filtfilt(b_lp, a_lp, filtered_data_hp, padlen=3 * (max(len(b_lp), len(a_lp)) - 1), axis=-1)

    # Differentiation to enhance R-peaks
    diff_data = np.diff(filtered_data, prepend=0)

    # Squaring to further emphasize R-peaks
    squared_data = diff_data ** 2

    # Moving average integration
    window_length = round(params['window_length'] * fs)
    window = np.ones(window_length) / window_length
    integrated_data = convolve(squared_data, window, mode='same')

    # Amplitude thresholding
    threshold = params['threshold_ratio'] * np.max(integrated_data)

    # Refractory period to avoid detecting multiple peaks within a short duration
    refractory_half_period = round(params['refractory_period'] * fs / 2)

    # Find the local peaks that satisfy both amplitude and width condition
    peaks = np.zeros(data.shape[0])
    for k in range(1, data.shape[0] + 1):
        search_win_start = max(1, k - refractory_half_period)
        search_win_stop = min(data.shape[0], k + refractory_half_period + 1)
        if threshold < integrated_data[k - 1] == np.max(integrated_data[search_win_start:search_win_stop]):
            peaks[k - 2] = 1
    # Get peak indexes
    peak_indexes = np.nonzero(peaks)[0] + 1

    return peaks, peak_indexes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_detection_pan_tompkins - R-peak detector based on Pan-Tompkins method.

    Args:
        data (numpy.ndarray): Vector of input ECG data
        fs (float): Sampling rate in Hz
        fc_low (float, optional): BP filter lower cutoff frequency in Hz (default: 5.0 Hz)
        fc_high (float, optional): BP filter upper cutoff frequency in Hz (default: 15.0 Hz)
        window_length (float, optional): Integration window length in seconds (default: 0.150 s)
        threshold_ratio (float, optional): Threshold ratio for peak detection (default: 0.2)
        refractory_period (float, optional): Refractory period in seconds (default: 0.2 s)

    Returns:
        peaks (numpy.ndarray): Vector of R-peak impulse train
        peak_indexes (numpy.ndarray): Vector of R-peak indexes

    References:
        Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
        Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532

    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_pan_tompkins.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    args = parser.parse_args()
