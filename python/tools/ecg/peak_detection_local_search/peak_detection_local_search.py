import sys

import numpy as np

sys.path.append('../peak_detection_simple')
from peak_detection_simple import peak_detection_simple


def peak_detection_local_search(x, ff, *args):
    """
    peak_detection_local_search - R-peak detector based on local max/min search
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
        July 2023: Translated to Python from Matlab (peak_detection_local_search.m)

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    # Set default values for optional arguments
    if len(args) > 0 and args[0] is not None:
        flag = args[0]
    else:
        flag = (np.abs(np.max(x))) > np.abs(np.min(x))

    if len(args) > 1 and args[1] is not None:
        num_rounds = args[1]
        if num_rounds > 3:
            raise ValueError('Number of R-peak detections not supported')
    else:
        num_rounds = 1

    if len(args) > 2 and args[2] is not None:
        hr_update_fraction = args[2]
    else:
        hr_update_fraction = 1.05

    if len(args) > 3 and args[3] is not None:
        omit_close_peaks = args[3]
    else:
        omit_close_peaks = 0

    # Perform peak detection
    peaks, peak_indexes = peak_detection_simple(x, ff, flag)

    # Perform additional iterations if specified
    if num_rounds > 1:
        for k in range(1, num_rounds):
            rr_intervals = np.diff(peak_indexes)
            ff = hr_update_fraction / np.median(rr_intervals)  # refined heart rate (in Hz) used for R-peak detection
            peaks, peak_indexes = peak_detection_simple(x, ff, flag, omit_close_peaks)
    return peaks, peak_indexes
