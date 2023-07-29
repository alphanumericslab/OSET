import argparse

import numpy as np


def peak_detection_simple(x, ff, flag=0, omit_close_peaks=0):
    """
    peak_detection_simple - Internal R-peak detector function
    Syntax: peaks, peak_indexes = peak_detection_simple(x, ff, flag, omit_close_peaks)

      Inputs:
          x: Vector of input data
          ff: Approximate ECG beat-rate in Hertz
          flag: Search for positive (flag=1) or negative (flag=0) peaks
          omit_close_peaks: omit close peaks after main peak detection
          (true/1) or not(false/0). Default is 0

      Outputs:
          peaks: Vector of R-peak impulse train
          peak_indexes: Vector of R-peak indexes

    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_simple.m)

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    n = len(x)
    peaks = np.zeros(n)
    rng = int(np.floor(0.5 / ff))
    if flag:
        for j in range(n):
            # Determine the index range for peak search
            if rng < j < n - rng:
                index = slice(j - rng - 1, j + rng + 1)
            elif j > rng:
                index = slice((n - 2 * rng) - 1, n + 1)
            else:
                index = slice(0, 2 * rng + 1)

            if np.max(x[index]) == x[j]:
                peaks[j] = 1
    else:
        for j in range(n):
            # Determine the index range for peak search
            if rng < j < n - rng:
                index = slice(j - rng - 1, j + rng + 1)
            elif j > rng:
                index = slice(n - 2 * rng - 1, n + 1)
            else:
                index = slice(0, 2 * rng + 1)

            if np.min(x[index]) == x[j]:
                peaks[j] = 1

    if omit_close_peaks:
        I = np.where(peaks)[0]
        d = np.diff(I)
        peaks[(I[np.where(d < rng)[0]])] = 0
    peak_indexes = np.where(peaks)[0] + 1
    return peaks, peak_indexes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_detection_simple - Internal R-peak detector function
    Syntax: peaks, peak_indexes = peak_detection_simple(x, ff, flag, omit_close_peaks)

      Inputs:
          x: Vector of input data
          ff: Approximate ECG beat-rate in Hertz
          flag: Search for positive (flag=1) or negative (flag=0) peaks
          omit_close_peaks: omit close peaks after main peak detection
          (true/1) or not(false/0). Default is 0

      Outputs:
          peaks: Vector of R-peak impulse train
          peak_indexes: Vector of R-peak indexes

    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_simple.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    args = parser.parse_args()
