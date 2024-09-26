import argparse

import numpy as np
from oset.ecg.peak_detection.peak_det_local_search import (
    peak_det_local_search,
)
from scipy.signal import lfilter


def peak_det_matched_filter_robust(ref, fs, h, fmax, itr=1):
    """
    peak_det_matched_filter_robust: R-peak detector based on a matched filter

    Syntax:
      peaks, r = peak_det_matched_filter_robust(ref, fs, h, fmax, itr)

    Args:
      ref (numpy.ndarray): vector of input data.
      fs (int): sampling rate.
      h (numpy.ndarray): template waveform.
      fmax (float): maximum expected frequency of the R-peaks.
      itr (int, optional): number of iterations to run the post-matched filter peak detector. Default is 1.


    Returns:
        peaks:  Vector of R-peak impulse train.
        r:      Filtered output after matched filtering.

      Revision History:
      July 2023: Translated to Python from Matlab (peak_det_matched_filter_robust.m)
      June 2024: Reworked the code to match the updated Matlab version

        Amulya, 2023
        The Open-Source Electrophysiological Toolbox
        https://github.com/alphanumericslab/OSET
    """
    N = len(ref)
    L = len(h)

    h = h[::-1]  # Reverse the template waveform

    w = L // 2

    # Matched filtering
    r = lfilter(h, 1, np.concatenate((ref, np.zeros(w - 1))))

    # Trim the filtered output to match the length of ref
    r = r[w - 1 : N + w - 1]
    peaks = peak_det_local_search(r, fmax / fs, 1, itr)[0]
    return peaks, r


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_det_matched_filter_robust: R-peak detector based on a matched filter

    Syntax:
      peaks, r = peak_det_matched_filter_robust(ref, fs, h, fmax, itr)

    Args:
      ref (numpy.ndarray): vector of input data.
      fs (int): sampling rate.
      h (numpy.ndarray): template waveform.
      fmax (float): maximum expected frequency of the R-peaks.
      itr (int, optional): number of iterations to run the post-matched filter peak detector. Default is 1.


    Returns:
        peaks:  Vector of R-peak impulse train.
        r:      Filtered output after matched filtering.

      Revision History:
      July 2023: Translated to Python from Matlab (peak_det_matched_filter_robust.m)
      June 2024: Reworked the code to match the updated Matlab version

        Amulya, 2023
        The Open-Source Electrophysiological Toolbox
        https://github.com/alphanumericslab/OSET
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
