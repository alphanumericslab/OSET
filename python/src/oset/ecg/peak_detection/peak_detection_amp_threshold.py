import argparse

import numpy as np
from oset.ecg.peak_detection.peak_detection_simple import peak_detection_simple


def peak_detection_amp_threshold(x, ff, th, *args):
    """
    peak_detection_amp_threshold - R-peak detector based on max search and level thresholding

    Syntax: peaks, peak_indexes = peak_detection_amp_threshold(x, ff, th, *args)

    Inputs:
      x: Vector of input data.
      ff: Approximate ECG beat-rate in Hertz, normalized by the sampling frequency.
      th: Peaks smaller than this fraction of the max peak amplitude are neglected.
      flag (optional): Search for positive (flag=1) or negative (flag=0) peaks. By default,
                  the maximum absolute value of the signal determines the peak sign.

    Outputs:
      peaks: Vector of R-peak impulse train.
      peak_indexes: Indexes of the detected R-peaks.

    Notes:
    - The R-peaks are found from a peak search in windows of length N, where
      N corresponds to the R-peak period calculated from the given ff. R-peaks
      with periods smaller than N/2 or greater than N are not detected.
    - The signal baseline wander is recommended to be removed before the
      R-peak detection.

    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_amp_threshold.m)

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    omit_close_peaks = 0
    if len(args) == 1:
        flag = args[0]
    else:
        flag = (np.abs(np.max(x))) > np.abs(np.min(x))
    peaks = peak_detection_simple(x, ff, flag, omit_close_peaks)[0]

    I = np.nonzero(peaks)[0]
    mx = np.max(np.abs(x[I]))
    J = np.abs(x[I]) < th * mx
    peaks[I[J]] = 0
    peak_indexes = np.where(peaks)[0] + 1
    return peaks, peak_indexes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_detection_amp_threshold - R-peak detector based on max search and level thresholding

    Syntax: peaks, peak_indexes = peak_detection_amp_threshold(x, ff, th, *args)

    Inputs:
      x: Vector of input data.
      ff: Approximate ECG beat-rate in Hertz, normalized by the sampling frequency.
      th: Peaks smaller than this fraction of the max peak amplitude are neglected.
      flag (optional): Search for positive (flag=1) or negative (flag=0) peaks. By default,
                  the maximum absolute value of the signal determines the peak sign.

    Outputs:
      peaks: Vector of R-peak impulse train.
      peak_indexes: Indexes of the detected R-peaks.

    Notes:
    - The R-peaks are found from a peak search in windows of length N, where
      N corresponds to the R-peak period calculated from the given ff. R-peaks
      with periods smaller than N/2 or greater than N are not detected.
    - The signal baseline wander is recommended to be removed before the
      R-peak detection.
    
    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_amp_threshold.m)
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
