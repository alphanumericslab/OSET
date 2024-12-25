import argparse

import numpy as np
from oset.ecg.peak_detection.peak_det_simple import peak_det_simple


def peak_det_local_search(
    x, ff, flag=None, num_rounds=1, hr_update_fraction=1.05, omit_close_peaks=0
):
    """
    peak_det_local_search - R-peak detector based on local max/min search

    Syntax:
      peaks, peak_indexes = peak_det_local_search(x, ff, *args)

    Args:
      x (numpy.ndarray): vector of input data.
      ff (float): approximate ECG beat-rate in Hertz, normalized by the sampling frequency.
      flag (int, optional): search for positive (flag=1) or negative (flag=0) peaks. By default,
                  the maximum absolute value of the signal determines the peak sign.
      num_rounds (int, optional): number of iterations to find the R-peaks, up to 3
                  (each time updating the expected R-peak rates). Default is 1 (no iterations).
      hr_update_fraction (float, optional): fraction to update the heart rate. Default is 1.05.
      omit_close_peaks (int, optional): whether to omit close peaks. Default is 0.


    output:
    peaks: vector of R-peak impulse train
    peak_indexes: vector of R-peak indexes

    Notes:
    - The R-peaks are found from a peak search in windows of length N; where
    N corresponds to the R-peak period calculated from the given f. R-peaks
    with periods smaller than N/2 or greater than N are not detected.
    - The signal baseline wander is recommended to be removed before the
    R-peak detection

    Revision History:
        July 2023: Translated to Python from Matlab (peak_det_local_search.m)
        July 2024: Update arguments to be more pythonic.

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    # Set default values for optional arguments
    if flag is None:
        flag = (np.abs(np.max(x))) > np.abs(np.min(x))

    if num_rounds > 3:
        raise ValueError("Number of R-peak detections not supported")

    # Perform peak detection
    peaks, peak_indexes = peak_det_simple(x, ff, flag, omit_close_peaks)

    # Perform additional iterations if specified
    if num_rounds > 1:
        for k in range(1, num_rounds):
            rr_intervals = np.diff(peak_indexes)
            ff = hr_update_fraction / np.median(
                rr_intervals
            )  # refined heart rate (in Hz) used for R-peak detection
            peaks, peak_indexes = peak_det_simple(x, ff, flag, omit_close_peaks)
    return peaks, peak_indexes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_det_local_search - R-peak detector based on local max/min search
    peaks, peak_indexes = PeakDetection(x,f,flag, num_rounds),

    inputs:
    x: vector of input data
    f: approximate ECG beat-rate in Hertz, normalized by the sampling frequency
    flag: search for positive (flag=1) or negative (flag=0) peaks. By default
    the maximum absolute value of the signal, determines the peak sign.
    num_rounds: the number of iterations to find the R-peaks, up to 3
    (everytime updating the expected R-peak rates). Default = 1 (no iterations)

    output:
    peaks: vector of R-peak impulse train
    peak_indexes: vector of R-peak indexes

    Notes:
    - The R-peaks are found from a peak search in windows of length N; where
    N corresponds to the R-peak period calculated from the given f. R-peaks
    with periods smaller than N/2 or greater than N are not detected.
    - The signal baseline wander is recommended to be removed before the
    R-peak detection

    Revision History:
        July 2023: Translated to Python from Matlab (peak_det_local_search.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
