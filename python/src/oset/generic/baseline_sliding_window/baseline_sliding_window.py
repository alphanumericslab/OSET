import argparse
import numpy as np


def baseline_sliding_window(x, l, approach):
    """
    baseline_sliding_window - Baseline wander extraction from biomedical recordings using a single stage of median or moving average filtering.

    Syntax: b = baseline_sliding_window(x, L, approach)

    Inputs:
      x: Vector or matrix of noisy data (channels x samples).
      L: Averaging window length (in samples).
      approach: Approach for baseline extraction.
          - 'md': Median filtering.
          - 'mn': Moving average.

    Output:
      b: Vector or matrix of baseline wanders (channels x samples).

    Revision History:
      Aug 2023 : Translated to Python from Matlab(baseline_sliding_window.m)

    Amulya, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    x = np.array(x, dtype=np.double)
    if x.ndim == 1:
        x = np.array([x], dtype=np.double)
    n = x.shape[1]
    b = np.zeros_like(x)
    flen = int(np.floor(l / 2))

    if approach == "mn":
        # moving average filter
        for j in range(1, n + 1):
            index = slice(max(0, j - flen - 1), min(n, j + flen))
            b[:, j - 1] = np.mean(x[:, index], axis=1)
    elif approach == "md":
        # median filter
        for j in range(1, n + 1):
            index = slice(max(0, j - flen - 1), min(n, j + flen))
            b[:, j - 1] = np.median(x[:, index], axis=1)
    return b


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    baseline_sliding_window - Baseline wander extraction from biomedical recordings using a single stage of median or moving average filtering.

    Syntax: b = baseline_sliding_window(x, L, approach)

    Inputs:
      x: Vector or matrix of noisy data (channels x samples).
      L: Averaging window length (in samples).
      approach: Approach for baseline extraction.
          - 'md': Median filtering.
          - 'mn': Moving average.

    Output:
      b: Vector or matrix of baseline wanders (channels x samples).

    Revision History:
      Aug 2023 : Translated to Python from Matlab(baseline_sliding_window.m)

        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
