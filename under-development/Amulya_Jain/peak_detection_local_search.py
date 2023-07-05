import numpy as np


def peak_detection(x, ff, *args):
    """
    peaks = PeakDetection(x,f,flag, num_rounds),
    R-peak detector based on max search

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
    """

    def peak_detection_internal(x, ff, flag):
        n = len(x)
        peaks = np.zeros(n)

        th = 0.5
        rng = int(np.floor(th / ff))

        if flag:
            for j in range(n):
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
                if rng < j < n - rng:
                    index = slice(j - rng - 1, j + rng + 1)
                elif j > rng:
                    index = slice(n - 2 * rng - 1, n + 1)
                else:
                    index = slice(0, 2 * rng + 1)

                if np.min(x[index]) == x[j]:
                    peaks[j] = 1

        # remove fake peaks
        I = np.where(peaks)[0]
        d = np.diff(I)
        peaks[(I[np.where(d < rng)[0]])] = 0
        peak_indexes = np.where(peaks)[0] + 1

        return peaks, peak_indexes

    # _____________________________________________________________________________________#
    if len(args) > 2 and args[0] is not None:
        flag = args[0]
    else:
        flag = (np.abs(np.max(x))) > np.abs(np.min(x))

    if len(args) > 3 and args[1] is not None:
        num_rounds = args[1]
        if num_rounds > 3:
            raise ValueError('Number of R-peak detections not supported')
    else:
        num_rounds = 1
    peaks, peak_indexes = peak_detection_internal(x, ff, flag)
    if num_rounds > 1:
        for k in range(1, num_rounds):
            rr_intervals = np.diff(peak_indexes)
            ff = 1.05 / np.median(rr_intervals)  # refined heart rate (in Hz) used for R-peak detection
            peaks, peak_indexes = peak_detection_internal(x, ff, flag)
    return peaks, peak_indexes
