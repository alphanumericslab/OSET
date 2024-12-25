import argparse

import numpy as np


def peak_det_simple(x, ff, flag, omit_close_peaks=0):
    """
    peak_det_simple - Internal R-peak detector function
    Syntax: peaks, peak_indexes = peak_det_simple(x, ff, flag, omit_close_peaks)

    Args:
        x (numpy.ndarray): vector of input data (ECG signal).
        ff (float): approximate ECG beat-rate in Hertz (normalized by the sampling frequency).
        flag (int): Search for positive (flag=1) or negative (flag=0) peaks.
        omit_close_peaks (int, optional): If non-zero (True), omit close peaks after main peak detection.
            Default is 0 (False).

      Returns:
          peaks: Vector of R-peak impulse train
          peak_indexes: Vector of R-peak indexes

    Revision History:
        July 2023: Translated to Python from Matlab (peak_det_simple.m)
        June 2024: Reworked the code to match the updated Matlab version

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    n = len(x)

    rng = int(np.floor(0.5 / ff))

    peaks = np.zeros(n)

    for j in range(1, n + 1):
        index = None
        # Determine the index range for peak search
        if rng < j < n - rng:
            index = slice(j - rng - 1, j + rng)
        elif j > rng:
            index = slice(j - rng - 1, n)
        else:
            index = slice(0, 2 * rng)
        if max(x[index]) == x[j - 1] and np.abs(x[index]).sum() > 0:
            peaks[j - 1] = 1

    nadirs = np.zeros(n)

    for j in range(1, n + 1):
        index = None
        # Determine the index range for peak search
        if rng < j < n - rng:
            index = slice(j - rng - 1, j + rng)
        elif j > rng:
            index = slice(n - 2 * rng - 1, n)
        else:
            index = slice(0, 2 * rng)
        if min(x[index]) == x[j - 1] and np.abs(x[index]).sum() > 0:
            nadirs[j - 1] = 1

    if flag == 0:
        peaks = nadirs
    elif flag == 2:
        if np.abs(np.median(x[peaks == 1])) < np.abs(np.median(x[nadirs == 1])):
            peaks = nadirs

    if omit_close_peaks:
        I = np.where(peaks)[0]
        d = np.diff(I)
        peaks[(I[d < rng])] = 0
    peak_indexes = np.where(peaks)[0] + 1
    return peaks, peak_indexes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_det_simple - Internal R-peak detector function
    Syntax: peaks, peak_indexes = peak_det_simple(x, ff, flag, omit_close_peaks)

    Args:
        x (numpy.ndarray): vector of input data (ECG signal).
        ff (float): approximate ECG beat-rate in Hertz (normalized by the sampling frequency).
        flag (int): Search for positive (flag=1) or negative (flag=0) peaks.
        omit_close_peaks (int, optional): If non-zero (True), omit close peaks after main peak detection.
            Default is 0 (False).

      Returns:
          peaks: Vector of R-peak impulse train
          peak_indexes: Vector of R-peak indexes

    Revision History:
        July 2023: Translated to Python from Matlab (peak_det_simple.m)
        June 2024: Reworked the code to match the updated Matlab version

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
