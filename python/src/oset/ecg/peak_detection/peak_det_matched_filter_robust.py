import argparse

import numpy as np
from oset.ecg.peak_detection.peak_det_local_search import (
    peak_det_local_search,
)
from scipy.signal import lfilter


def peak_det_matched_filter_robust(ref, fs, h, fmax, itr=1):
    """
    peak_det_matched_filter_robust: R-peak detector based on a matched filter

    Syntax: [peaks, r] = peak_det_matched_filter_robust(ref, fs, h, fmax, itr)

    Args:
        ref:    Vector of input data.s
        fs:     Sampling rate.
        h:      Template waveform.
        fmax:   Maximum expected frequency of the R-peaks.
        itr:    Number of iterations to run the post-matched filter peak detector (Default: 1, run peak detector only once)

    Returns:
        peaks:  Vector of R-peak impulse train.
        r:      Filtered output after matched filtering.

      Revision History:
      July 2023: Translated to Python from Matlab (peak_det_matched_filter_robust.m)

        Amulya, 2023
        The Open-Source Electrophysiological Toolbox
        https://github.com/alphanumericslab/OSET
    """
    n = len(ref)
    length = len(h)

    h = np.flip(h)  # Reverse the template waveform
    w = int(np.floor(length / 2))
    r = lfilter(h, 1, np.concatenate((ref, np.zeros(w - 1))))
    r = r[w - 1 : w + n]
    peaks = peak_det_local_search(r, fmax / fs, 1, itr)[0]
    return peaks, r


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        peak_det_matched_filter_robust: R-peak detector based on a matched filter

        Syntax: peaks, r = peak_det_matched_filter_robust(ref, fs, h, fmax, itr)

        Inputs:
          ref:    Vector of input data.
          fs:     Sampling rate.
          h:      Template waveform.
          fmax:   Maximum expected frequency of the R-peaks.
          itr:    Number of iterations to run the post-matched filter peak detector (Default: 1, run peak detector only once)

        Outputs:
          peaks:  Vector of R-peak impulse train.
          r:      Filtered output after matched filtering.

          Revision History:
          July 2023: Translated to Python from Matlab (peak_det_matched_filter_robust.m)

            Amulya, 2023
            The Open-Source Electrophysiological Toolbox
            https://github.com/alphanumericslab/OSET
        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
