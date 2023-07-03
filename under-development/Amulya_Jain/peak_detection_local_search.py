import numpy as np


def peak_detection(x, ff, *args):
    """
 peak_detection_local_search - R-peak detector based on local max/min search

   [peaks, peak_indexes] = peak_detection_local_search(x, f, flag, num_rounds, hr_update_fraction, omit_close_peaks)

   Inputs:
       x: Vector of input data
       f: Approximate ECG beat-rate in Hertz, normalized by the sampling frequency
       flag: Optional. Search for positive (flag=1) or negative (flag=0) peaks.
             By default, the maximum absolute value of the signal determines the peak sign.
       num_rounds: Optional. Number of iterations to find the R-peaks, up to 3.
                   Default is 1 (run peak detection only once).
       hr_update_fraction: Optional. median HR is multiplied by this
                   fraction in each iteration, if num_rounds > 1. Default
                   is 1.05
       omit_close_peaks: omit 'too-close' peaks after main peak detection (true/1) or not(false/0) (Default is false)

   Outputs:
       peaks: Vector of R-peak impulse train
       peak_indexes: Vector of R-peak indexes

   Notes:
       - The R-peaks are found from a peak search in windows of length N, where
         N corresponds to the R-peak period calculated from the given f. R-peaks
         with periods smaller than N/2 or greater than N are not detected.
       - It is recommended to remove the signal baseline wander before R-peak detection.
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