import argparse
from math import floor

import numpy as np
from oset.generic.baseline_sliding_window.baseline_sliding_window import (
    baseline_sliding_window,
)


def peak_det_adaptive_hr(
    x, ff, fs, flag=None, th=0.5, th2=0.8, rejection_threshold=0.3
):
    """
    peak_detection_adaptive_hr - R-peak detector based on max search over sliding window with adaptive width
    Syntax:
      peaks, peak_indexes = peak_detection_adaptive_hr(x, ff, fs, flag, [th, th2], rejection_threshold)

    Method:
      Detects R-peaks in the ECG signal x using a max search algorithm over a sliding window
      with adaptive width. The function returns the binary impulse train of detected peaks
      (PEAKS) and their corresponding indexes (PEAK_INDEXES) in the input signal x.

    Inputs:
      x (numpy.ndarray): vector of the input ECG signal.
      ff (float): approximate ECG beat-rate in Hz.
      fs (float): sampling frequency in Hz.
      flag (int, optional): search for positive (flag=1) or negative (flag=0) peaks.
          By default, the maximum absolute value of the signal determines the peak sign.
      th (float, optional): threshold used for the first round of peak detection (default: 0.5).
      th2 (float, optional): threshold used for the second round of peak detection (default: 0.8).
      rejection_threshold (float, optional): threshold for rejecting fake peaks (default: 0.3).

    Outputs:
      peaks: Binary impulse train indicating the detected R-peaks.
      peak_indexes: Indexes of the detected R-peaks in the input signal.

    Note:
      - The signal baseline wander should be removed before R-peak detection.

    Revision History:
      June 2024: Translated to Python from Matlab(peak_detection_adaptive_hr.m)

    Amulya Jain, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    if flag is None:
        flag = abs(max(x)) > abs(min(x))

    n = len(x)
    peaks0 = np.zeros(n)
    peaks = np.zeros(n)

    rng = floor(th / (ff / fs))

    # First round peak detection
    if flag:
        for j in range(1, n + 1):
            # Determine the index range for peak search
            if rng < j < n - rng:
                index = slice(j - rng - 1, j + rng)
            elif j > rng:
                index = slice((n - 2 * rng) - 1, n)
            else:
                index = slice(0, 2 * rng)

            if np.max(x[index]) == x[j - 1]:
                peaks0[j - 1] = 1
    else:
        for j in range(1, n + 1):
            # Determine the index range for peak search
            if rng < j < n - rng:
                index = slice(j - rng - 1, j + rng)
            elif j > rng:
                index = slice(n - 2 * rng - 1, n)
            else:
                index = slice(0, 2 * rng)

            if np.min(x[index]) == x[j - 1]:
                peaks0[j - 1] = 1

    # Smooth the RR intervals
    ii = np.array(np.where(peaks0), dtype=np.double) + 1
    rr_intervals = np.diff(ii)
    rr_intervals_smoothed = baseline_sliding_window(rr_intervals, 3, "md")
    rr_intervals_smoothed2 = baseline_sliding_window(rr_intervals_smoothed, 3, "mn")
    ff = fs / np.concatenate(
        (np.array([rr_intervals_smoothed2[0]]), rr_intervals_smoothed2)
    )

    # Handle edge cases
    if (ii[0] > 1).all():
        ii = np.insert(ii, 0, 1)
        ff = np.insert(ff, 0, np.mean(ff[0:2]))
    if (ii[-1] < n).all():
        ii = np.append(ii, n)
        ff = np.append(ff, np.mean(ff[-2:]))

    # Interpolate the beat-rate for all samples
    ff_interpolated = np.interp(np.arange(1, n + 1), ii, ff)
    rng2 = np.floor(th2 / (ff_interpolated / fs))  # Small difference

    # Second round peak detection
    if flag:
        for j in range(1, n + 1):
            # Determine the index range for peak search
            if rng2[j - 1] < j < n - rng2[j - 1]:
                index = slice(int(j - rng2[j - 1] - 1), int(j + rng2[j - 1]))
            elif j > rng2[j - 1]:
                index = slice(int((n - 2 * rng2[j - 1]) - 1), int(n))
            else:
                index = slice(0, int(2 * rng2[j - 1]))

            if np.max(x[index]) == x[j - 1]:
                peaks[j - 1] = 1
    else:
        for j in range(1, n + 1):
            # Determine the index range for peak search
            if rng2[j - 1] < j < n - rng2[j - 1]:
                index = slice(int(j - rng2[j - 1] - 1), int(j + rng2[j - 1]))
            elif j > rng2[j - 1]:
                index = slice(int(n - 2 * rng2[j - 1] - 1), n)
            else:
                index = slice(0, int(2 * rng2[j - 1]))

            if np.min(x[index]) == x[j - 1]:
                peaks[j - 1] = 1

    # Remove fake peaks
    i = np.array(np.where(peaks), dtype=np.int64) + 1
    i = i[0]
    peak_amps = np.median(x[i - 1])
    j = np.array(
        np.abs(x[i - 1]) < rejection_threshold * np.abs(peak_amps), dtype=np.bool_
    )
    peaks[i[j] - 1] = 0

    peak_indexes = np.array(np.where(peaks), dtype=np.int64) + 1
    peak_indexes = peak_indexes[0]
    return peaks, peak_indexes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_detection_adaptive_hr - R-peak detector based on max search over sliding window with adaptive width
    Syntax:
      peaks, peak_indexes = peak_detection_adaptive_hr(x, ff, fs, flag, [th, th2], rejection_threshold)

    Method:
      Detects R-peaks in the ECG signal x using a max search algorithm over a sliding window
      with adaptive width. The function returns the binary impulse train of detected peaks
      (PEAKS) and their corresponding indexes (PEAK_INDEXES) in the input signal x.

    Inputs:
      x (numpy.ndarray): vector of the input ECG signal.
      ff (float): approximate ECG beat-rate in Hz.
      fs (float): sampling frequency in Hz.
      flag (int, optional): search for positive (flag=1) or negative (flag=0) peaks.
          By default, the maximum absolute value of the signal determines the peak sign.
      th (float, optional): threshold used for the first round of peak detection (default: 0.5).
      th2 (float, optional): threshold used for the second round of peak detection (default: 0.8).
      rejection_threshold (float, optional): threshold for rejecting fake peaks (default: 0.3).

    Outputs:
      peaks: Binary impulse train indicating the detected R-peaks.
      peak_indexes: Indexes of the detected R-peaks in the input signal.

    Note:
      - The signal baseline wander should be removed before R-peak detection.

    Revision History:
      June 2024: Translated to Python from Matlab(peak_detection_adaptive_hr.m)

    Amulya Jain, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
